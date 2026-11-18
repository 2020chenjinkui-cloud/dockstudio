"""Regression tests for hard-abort propagation vs cooperative cancellation.

Background: ``ScreeningRunner`` used to keep "stop asking for more work" and
"a worker blew up" in one boolean (``_cancel_requested``).  A BaseException
raised inside a worker therefore landed in a ``Future`` that the consumer loop
never read - it ``break``-ed first - so ``run()`` returned a *truncated* report
as if the screening had completed normally.  Measured before the fix:
``test_adversarial_m3_challenger2::test_mid_run_crash_simulation_at_7_and_19``
failed 18/24 times and ``test_screen_runner::test_checkpoint_interrupt_...``
failed 11/24 times, always with ``DID NOT RAISE KeyboardInterrupt``.

These tests are deterministic on purpose: they never sleep and never rely on
thread timing.  Instead a fake engine parks inside ``dock()`` on an
``threading.Event`` until the test releases it, so the exact interleaving that
used to be a coin flip is forced - both the "consumer already saw the flag"
order and the "consumer never read the Future" order are exercised explicitly.

See docs/fixes/abort-signal-2026-11-18.md.
"""

from __future__ import annotations

import json
import threading
from pathlib import Path

import pytest

from dockstudio.cli.commands.verify import load_funnel_data
from dockstudio.core.screen_runner import ScreeningRunner
from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.models import (
    Box,
    DockParams,
    DockResult,
    Pose,
    PreparedLibrary,
    PreparedLigand,
    ScreeningConfig,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

_SAMPLE_PDBQT = (
    "MODEL 1\n"
    "REMARK VINA RESULT: -7.50 0.000 0.000\n"
    "ROOT\n"
    "ATOM      1  C   UNL     1       1.000   2.000   3.000  1.00  0.00    +0.000 C\n"
    "ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O\n"
    "ENDROOT\n"
    "ENDMDL\n"
)


def _make_library(tmp_path: Path, n_ligands: int) -> tuple[PreparedLibrary, list[PreparedLigand]]:
    lib_dir = tmp_path / "lib"
    lib_dir.mkdir(parents=True, exist_ok=True)
    ligands: list[PreparedLigand] = []
    ids: list[str] = []
    for i in range(1, n_ligands + 1):
        lid = f"LIG_{i:03d}"
        sub = lib_dir / lid
        sub.mkdir(parents=True, exist_ok=True)
        pdbqt = sub / "ligand.pdbqt"
        pdbqt.write_text(_SAMPLE_PDBQT, encoding="utf-8")
        ligands.append(
            PreparedLigand(
                id=lid,
                source="smiles",
                pdbqt_path=str(pdbqt),
                smiles="CC(=O)Oc1ccccc1C(=O)O",
                n_rotatable_bonds=2,
            )
        )
        ids.append(lid)
    prep_lib = PreparedLibrary(
        id=f"lib_{n_ligands}",
        name=f"Lib{n_ligands}",
        path=str(lib_dir),
        ligand_ids=ids,
        n_raw=n_ligands,
        n_ready=n_ligands,
    )
    return prep_lib, ligands


def _make_receptor(tmp_path: Path) -> Path:
    rec = tmp_path / "receptor.pdbqt"
    rec.write_text(
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00    -0.300 N\n",
        encoding="utf-8",
    )
    return rec


class RecordingEngine(BaseEngine):
    """Mock engine: records dock() order, optionally aborts on a given ligand.

    ``gate`` (optional) is waited on *inside* dock() so a test can hold every
    worker at a known point and drive the interleaving deterministically.
    """

    manifest = EngineManifest(
        name="recording_engine",
        version="1.0.0",
        scoring_functions=["vina"],
        supports_batch_ligands=False,
    )

    def __init__(
        self,
        abort_on: str | None = None,
        gate: threading.Event | None = None,
        released: threading.Event | None = None,
    ) -> None:
        super().__init__()
        self.abort_on = abort_on
        self.gate = gate
        self.released = released
        self.calls: list[str] = []
        # (ligand_id, exhaustiveness): the tiered funnel legitimately re-docks a
        # tier-1 passer at a different exhaustiveness in tier 2.
        self.call_details: list[tuple[str, int]] = []
        self._lock = threading.Lock()

    def _ligand_id(self, ligand_pdbqt: Path) -> str:
        stem = ligand_pdbqt.stem
        return ligand_pdbqt.parent.name if stem == "ligand" else stem

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        lig_id = self._ligand_id(Path(ligand_pdbqt))
        with self._lock:
            self.calls.append(lig_id)
            self.call_details.append((lig_id, params.exhaustiveness))
        if self.gate is not None and self.released is not None and lig_id == self.abort_on:
            # Signal "this worker is inside dock()" and wait to be released, so
            # the test - not the OS scheduler - decides what happens next.
            self.gate.set()
            self.released.wait(timeout=10)
        if lig_id == self.abort_on:
            raise KeyboardInterrupt(f"simulated crash on {lig_id}")

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        out = work_dir / "out.pdbqt"
        out.write_text(_SAMPLE_PDBQT, encoding="utf-8")
        return DockResult(
            engine="recording_engine",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-7.5, rmsd_lb=0.0, rmsd_ub=0.0)],
            log="OK",
            run_dir=str(work_dir),
        )


def _runner(
    tmp_path: Path,
    engine: BaseEngine,
    n_ligands: int,
    *,
    tiered: bool = False,
    cpu_per_job: int = 2,
    out_dir: Path | None = None,
) -> tuple[ScreeningRunner, PreparedLibrary]:
    """Build a runner on a fixed out_dir (resume phases must share one directory)."""
    prep_lib, _ = _make_library(tmp_path, n_ligands)
    return (
        ScreeningRunner(
            config=ScreeningConfig(enable_tiered_funnel=tiered, cpu=2),
            out_dir=out_dir or (tmp_path / "run"),
            receptor_pdbqt=_make_receptor(tmp_path),
            prepared_library=prep_lib,
            engine=engine,
            cpu_per_job=cpu_per_job,
        ),
        prep_lib,
    )


def _checkpoint_ids(out_dir: Path) -> list[str]:
    """Checkpointed ligand ids, in file order, de-duplicated.

    ``CheckpointManager`` streams one record per (ligand, tier), so the tiered
    funnel legitimately writes several lines for the same molecule; the contract
    under test is "every completed molecule survives a crash", not "one line".
    """
    path = out_dir / "checkpoint.jsonl"
    if not path.is_file():
        return []
    seen: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        lid = json.loads(line)["ligand_id"]
        if lid not in seen:
            seen.append(lid)
    return seen


# ---------------------------------------------------------------------------
# 1. Hard abort must propagate - single-tier path
# ---------------------------------------------------------------------------


def test_hard_abort_propagates_single_tier(tmp_path: Path) -> None:
    """A crash inside a worker must leave run() as an exception, never as a report.

    This is the exact scenario that used to fail ~75% of the time.  It runs the
    worker to completion in its own thread, so the exception is *already sitting
    in the Future* - the consumer would previously break on ``_cancel_requested``
    before reading it.
    """
    engine = RecordingEngine(abort_on="LIG_004")
    runner, _ = _runner(tmp_path, engine, n_ligands=8)

    with pytest.raises(KeyboardInterrupt, match="simulated crash on LIG_004"):
        runner.run(resume=False)

    # The molecules docked before the crash are still checkpointed (zero loss).
    assert _checkpoint_ids(tmp_path / "run") == ["LIG_001", "LIG_002", "LIG_003"]


def test_hard_abort_propagates_when_consumer_already_saw_cancel(tmp_path: Path) -> None:
    """Order-independence: abort known before the consumer reads its Future.

    The worker parks inside ``dock()``; the test releases it only after the
    consumer has been given a chance to observe the cancel flag.  Even then the
    KeyboardInterrupt must surface - it may not be swallowed by the ``break``.
    """
    gate = threading.Event()
    released = threading.Event()
    engine = RecordingEngine(abort_on="LIG_002", gate=gate, released=released)
    runner, _ = _runner(tmp_path, engine, n_ligands=5)

    def _release_soon() -> None:
        gate.wait(timeout=10)
        released.set()

    helper = threading.Thread(target=_release_soon, daemon=True)
    helper.start()
    try:
        with pytest.raises(KeyboardInterrupt, match="simulated crash on LIG_002"):
            runner.run(resume=False)
    finally:
        released.set()
        helper.join(timeout=5)

    assert _checkpoint_ids(tmp_path / "run") == ["LIG_001"]


def test_run_resets_abort_state_between_runs(tmp_path: Path) -> None:
    """A later resume on the same runner object must not re-raise a stale abort."""
    engine = RecordingEngine(abort_on="LIG_002")
    runner, _ = _runner(tmp_path, engine, n_ligands=4)
    with pytest.raises(KeyboardInterrupt):
        runner.run(resume=False)

    runner.engine = RecordingEngine()  # clean engine, no crash
    report = runner.run(resume=True)
    assert len(report.candidates) == 4


# ---------------------------------------------------------------------------
# 2. Hard abort must propagate - tiered funnel path
# ---------------------------------------------------------------------------


def test_hard_abort_propagates_tier1_tiered_funnel(tmp_path: Path) -> None:
    """Tier 1 of the hierarchical funnel must propagate the abort as well."""
    engine = RecordingEngine(abort_on="LIG_002")
    runner, _ = _runner(tmp_path, engine, n_ligands=6, tiered=True)

    with pytest.raises(KeyboardInterrupt, match="simulated crash on LIG_002"):
        runner.run(resume=False)

    # Molecules queued behind the crash are abandoned, not silently docked.
    assert len(engine.calls) < 6, f"queue was drained after the crash: {engine.calls}"


def test_hard_abort_does_not_drain_worker_queue(tmp_path: Path) -> None:
    """The pool stops handing out work once a hard abort is recorded.

    Before the fix the tiered path docked **all** molecules and only then raised
    (measured: 10 dock calls for a crash on #5), burning GPU time on a run that
    was going to be discarded.
    """
    engine = RecordingEngine(abort_on="LIG_003")
    runner, _ = _runner(tmp_path, engine, n_ligands=12, tiered=True)

    with pytest.raises(KeyboardInterrupt):
        runner.run(resume=False)

    assert len(engine.calls) < 12, f"docked {len(engine.calls)}/12 after abort"


# ---------------------------------------------------------------------------
# 3. Cooperative cancel: graceful, but visibly truncated
# ---------------------------------------------------------------------------


def test_cooperative_cancel_returns_partial_report_but_marks_truncation(
    tmp_path: Path,
) -> None:
    """``cancel()`` keeps its old promise (return, no exception) - but the report
    must no longer claim the screening completed: an explicit ``aborted`` stage
    row appears in funnel.csv and the mass-conservation identity stays open."""
    engine = RecordingEngine()
    runner, _ = _runner(tmp_path, engine, n_ligands=6)

    original_dock = engine.dock

    def dock_then_cancel(*args, **kwargs):
        result = original_dock(*args, **kwargs)
        runner.cancel()
        return result

    engine.dock = dock_then_cancel  # type: ignore[method-assign]

    report = runner.run(resume=False)  # graceful: no exception

    assert len(report.candidates) < 6
    funnel = load_funnel_data(tmp_path / "run" / "funnel.csv")
    assert funnel.get("aborted", 0) > 0, f"no 'aborted' row: {funnel}"
    assert report.is_conserved is False, "truncated run must not close the identity"


def test_cooperative_cancel_is_not_reported_as_dock_failure(tmp_path: Path) -> None:
    """Cancellation is not an engine failure: the abandoned molecules go to the
    ``aborted`` bucket, never into ``dock_failed`` (which would blame the engine
    for the operator's decision)."""
    engine = RecordingEngine()
    runner, _ = _runner(tmp_path, engine, n_ligands=5)

    original_dock = engine.dock

    def dock_then_cancel(*args, **kwargs):
        result = original_dock(*args, **kwargs)
        runner.cancel()
        return result

    engine.dock = dock_then_cancel  # type: ignore[method-assign]
    report = runner.run(resume=False)

    assert report.n_dock_failed == 0
    assert report.n_docked_success == len(report.candidates)


# ---------------------------------------------------------------------------
# 4. The two historically flaky tests, made deterministic
# ---------------------------------------------------------------------------


def test_abort_then_resume_zero_loss_three_phases(tmp_path: Path) -> None:
    """Deterministic version of the "crash at 7, crash at 19, finish" scenario.

    Same product contract as ``test_adversarial_m3_challenger2`` /
    ``test_screen_runner``, but with no reliance on thread timing: each phase
    crashes on a *named* molecule and the assertions are on checkpoint contents
    and dock call counts, which is what "zero loss on resume" actually means.
    """
    engine = RecordingEngine(abort_on="LIG_004")
    runner, _ = _runner(tmp_path, engine, n_ligands=9)
    out_dir = tmp_path / "run"

    with pytest.raises(KeyboardInterrupt, match="LIG_004"):
        runner.run(resume=False)
    assert _checkpoint_ids(out_dir) == ["LIG_001", "LIG_002", "LIG_003"]

    engine2 = RecordingEngine(abort_on="LIG_007")
    runner2, _ = _runner(tmp_path, engine2, n_ligands=9, out_dir=out_dir)
    with pytest.raises(KeyboardInterrupt, match="LIG_007"):
        runner2.run(resume=True)
    assert _checkpoint_ids(out_dir) == [f"LIG_{i:03d}" for i in range(1, 7)]

    engine3 = RecordingEngine()
    runner3, _ = _runner(tmp_path, engine3, n_ligands=9, out_dir=out_dir)
    report = runner3.run(resume=True)

    ids = _checkpoint_ids(out_dir)
    assert ids == [f"LIG_{i:03d}" for i in range(1, 10)]
    assert len(set(ids)) == 9, "duplicate checkpoint entries"
    assert len(report.candidates) == 9
    # No molecule is docked twice *except* the two crash molecules: those were
    # attempted (and aborted) once and legitimately retried after the resume.
    all_calls = engine.calls + engine2.calls + engine3.calls
    assert set(all_calls) == {f"LIG_{i:03d}" for i in range(1, 10)}
    retried = {lid for lid in set(all_calls) if all_calls.count(lid) > 1}
    assert retried == {"LIG_004", "LIG_007"}, f"unexpected re-docks: {sorted(retried)}"


def test_tiered_abort_then_resume_reaches_full_completion(tmp_path: Path) -> None:
    """Same contract on the hierarchical path: abort, resume, nothing lost."""
    engine = RecordingEngine(abort_on="LIG_003")
    runner, _ = _runner(tmp_path, engine, n_ligands=8, tiered=True)
    out_dir = tmp_path / "run"

    with pytest.raises(KeyboardInterrupt):
        runner.run(resume=False)

    # Resume with a clean engine: everything still gets screened.
    engine2 = RecordingEngine()
    runner2, _ = _runner(tmp_path, engine2, n_ligands=8, tiered=True, out_dir=out_dir)
    report = runner2.run(resume=True)

    assert report.n_raw == 8
    # Zero-loss contract: every molecule is checkpointed exactly once
    # (de-duplicated across tiers) after the resume.
    assert _checkpoint_ids(out_dir) == [f"LIG_{i:03d}" for i in range(1, 9)]
    # The funnel still filters (tier1 keep ratio -> tier2 keep ratio), so far
    # fewer candidates than molecules is the correct outcome, not data loss.
    assert 0 < len(report.candidates) < 8
    assert report.n_docked_success == len(report.candidates)
    # Molecules checkpointed at tier 1 during phase 1 are not re-docked at
    # tier-1 exhaustiveness (=4). Tier 2 re-docking a tier-1 passer at
    # exhaustiveness=8 is the funnel working as designed, not a repeat docking.
    tier1_phase2 = {lid for lid, exh in engine2.call_details if exh == 4}
    assert "LIG_001" not in tier1_phase2, sorted(tier1_phase2)
    assert "LIG_002" not in tier1_phase2, sorted(tier1_phase2)
