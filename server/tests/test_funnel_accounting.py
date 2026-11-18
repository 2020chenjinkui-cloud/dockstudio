"""Funnel accounting invariants: every molecule lands in exactly ONE bucket.

The invariant this module enforces is not "the numbers add up" but the stronger
statement the numbers are supposed to express:

    raw + aglycone_generated
      = tier1_filtered + tier2_filtered + dock_failed + aborted + docked_success

with the buckets **disjoint**.  Two opposite failures have been found in this
area, and both are covered here:

* t12-F2: a cooperative cancel left molecules that an earlier tier had already
  checkpointed, so they were absorbed by ``n_funnel_filtered`` and a truncated
  run closed (``is_conserved=True``) - "not conserved but claims to be".
* t15-F1: a *real* engine failure was counted twice (once in ``dock_failed``,
  again in ``tier1_filtered`` because a failed ligand has no usable affinity and
  therefore drops out of the tier's output count), so a normally completed run
  reported ``is_conserved=False`` and /verify T5 FAILED - "conserved but claims
  not to be".

Both are the same bug in opposite directions: a molecule was allowed to occupy
two buckets at once.
"""

from __future__ import annotations

import threading
from pathlib import Path

import pytest

from dockstudio.cli.commands.verify import _verify_t5, load_delivery_package, load_energy_table
from dockstudio.core.screen_runner import ScreeningRunner
from dockstudio.models import (
    Box,
    DockingError,
    DockParams,
    DockResult,
    Pose,
    PreparedLibrary,
    PreparedLigand,
    ScreeningConfig,
)

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


class FaultEngine:
    """Mock engine that can fail (DockingError) or crash (KeyboardInterrupt) per ligand."""

    manifest = type(
        "M",
        (),
        {
            "name": "fault_engine",
            "version": "1.0.0",
            "scoring_functions": ["vina"],
            "supports_batch_ligands": False,
        },
    )()

    def __init__(
        self,
        fail_on: set[str] | None = None,
        abort_on: str | None = None,
    ) -> None:
        self.fail_on = set(fail_on or ())
        self.abort_on = abort_on
        self.calls: list[str] = []
        self._lock = threading.Lock()

    def _ligand_id(self, ligand_pdbqt: Path) -> str:
        stem = Path(ligand_pdbqt).stem
        return Path(ligand_pdbqt).parent.name if stem == "ligand" else stem

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        lig_id = self._ligand_id(ligand_pdbqt)
        with self._lock:
            self.calls.append(lig_id)
        if self.abort_on is not None and lig_id == self.abort_on:
            raise KeyboardInterrupt(f"simulated crash on {lig_id}")
        if lig_id in self.fail_on:
            raise DockingError(f"simulated engine failure on {lig_id}")

        work_dir = Path(work_dir)
        work_dir.mkdir(parents=True, exist_ok=True)
        (work_dir / "out.pdbqt").write_text(_SAMPLE_PDBQT, encoding="utf-8")
        return DockResult(
            engine="fault_engine",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=-7.5, rmsd_lb=0.0, rmsd_ub=0.0)],
            log="OK",
            run_dir=str(work_dir),
        )


def _runner(
    tmp_path: Path,
    engine,
    n_ligands: int,
    *,
    keep: tuple[float, float] = (1.0, 1.0),
    out_dir: Path | None = None,
) -> ScreeningRunner:
    prep_lib, _ = _make_library(tmp_path, n_ligands)
    return ScreeningRunner(
        config=ScreeningConfig(
            enable_tiered_funnel=True,
            cpu=2,
            tier1_keep_ratio=keep[0],
            tier2_keep_ratio=keep[1],
        ),
        out_dir=out_dir or (tmp_path / "run"),
        receptor_pdbqt=_make_receptor(tmp_path),
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )


def _buckets(report, runner=None) -> dict[str, int]:
    """Scalar buckets plus the aborted count (which lives on the runner, not the
    report schema - see ScreeningRunner.abort_accounting)."""
    ids = runner.abort_accounting() if runner is not None else {"aborted": [], "failed": []}
    return {
        "tier1_filtered": report.n_tier1_filtered,
        "tier2_filtered": report.n_tier2_filtered,
        "funnel_filtered": report.n_funnel_filtered,
        "dock_failed": report.n_dock_failed,
        "aborted": len(ids["aborted"]),
        "docked_success": report.n_docked_success,
    }


# ---------------------------------------------------------------------------
# (a) real engine failures only -> the run DID complete -> must conserve
# ---------------------------------------------------------------------------


def test_tiered_single_engine_failure_conserves(tmp_path: Path) -> None:
    """One genuine DockingError must land in ``dock_failed`` and nowhere else.

    t15-F1: a failed ligand has no usable affinity, so it drops out of the tier's
    output count; deriving ``tier1_filtered`` from ``n_t1_in - n_t1_out`` counted
    it a second time and made a completed run report ``is_conserved=False``.
    """
    engine = FaultEngine(fail_on={"LIG_004"})
    runner = _runner(tmp_path, engine, 12)
    report = runner.run(resume=False)

    assert report.n_dock_failed == 1, _buckets(report, runner)
    assert report.n_funnel_filtered == 0, _buckets(report, runner)
    assert report.is_conserved is True, _buckets(report, runner)
    # The four buckets partition the pool exactly once.
    total = (
        report.n_funnel_filtered
        + report.n_dock_failed
        + len(runner.abort_accounting()["aborted"])
        + report.n_docked_success
    )
    assert total == report.n_docking_pool, _buckets(report, runner)


def test_tiered_two_engine_failures_conserve(tmp_path: Path) -> None:
    """Two failures must be counted twice, not four times."""
    engine = FaultEngine(fail_on={"LIG_003", "LIG_009"})
    runner = _runner(tmp_path, engine, 12)
    report = runner.run(resume=False)

    assert report.n_dock_failed == 2, _buckets(report, runner)
    assert report.n_funnel_filtered == 0, _buckets(report, runner)
    assert report.is_conserved is True, _buckets(report, runner)
    total = (
        report.n_funnel_filtered
        + report.n_dock_failed
        + len(runner.abort_accounting()["aborted"])
        + report.n_docked_success
    )
    assert total == report.n_docking_pool, _buckets(report, runner)


def test_tiered_failures_with_ratio_cutoff_conserve(tmp_path: Path) -> None:
    """Ratio cut-off + engine failures together: filtered must be the *ratio* only."""
    engine = FaultEngine(fail_on={"LIG_002"})
    runner = _runner(tmp_path, engine, 12, keep=(0.5, 1.0))
    report = runner.run(resume=False)

    # 12 in, 1 failed, 11 usable -> 50% keep = 6 pass. The other 5 are ratio losses.
    assert report.n_dock_failed == 1, _buckets(report, runner)
    assert report.n_tier1_filtered == 5, _buckets(report, runner)
    assert report.n_funnel_filtered == report.n_tier1_filtered + report.n_tier2_filtered
    assert report.is_conserved is True, _buckets(report, runner)
    total = (
        report.n_funnel_filtered
        + report.n_dock_failed
        + len(runner.abort_accounting()["aborted"])
        + report.n_docked_success
    )
    assert total == report.n_docking_pool, _buckets(report, runner)


def test_tiered_engine_failure_passes_verify_t5(tmp_path: Path) -> None:
    """End-to-end: /verify T5 must PASS for a completed run that had failures."""
    engine = FaultEngine(fail_on={"LIG_005"})
    runner = _runner(tmp_path, engine, 12)
    out_dir = tmp_path / "run"
    runner.run(resume=False)

    pkg = load_delivery_package(out_dir)
    energy_rows = load_energy_table(out_dir / "screening_results.csv")
    result = _verify_t5(pkg, energy_rows, None)

    assert result.passed is True, f"{result.status}: {result.observed}"
    assert "差值=0" in result.observed, result.observed


# ---------------------------------------------------------------------------
# (b) failure + cancel at the same time -> truncated, buckets still disjoint
# ---------------------------------------------------------------------------


def test_tiered_failure_and_cancel_keep_buckets_disjoint(tmp_path: Path) -> None:
    """A cancel must not swallow the failures, and vice versa.

    Invariants: the four buckets still partition the pool, and ``aborted`` and
    ``dock_failed`` never overlap (a molecule is either a genuine engine failure
    or it was never evaluated - not both).
    """
    prep_lib, _ = _make_library(tmp_path, 12)
    engine = FaultEngine(fail_on={"LIG_002", "LIG_003"})
    runner = ScreeningRunner(
        config=ScreeningConfig(enable_tiered_funnel=True, cpu=2),
        out_dir=tmp_path / "run",
        receptor_pdbqt=_make_receptor(tmp_path),
        prepared_library=prep_lib,
        engine=engine,
        cpu_per_job=2,
    )

    original_emit = ScreeningRunner._emit_progress
    state = {"seen": 0}

    def emit_then_cancel(self, tier, stage, current, total, lig_id, t_start):
        result = original_emit(self, tier, stage, current, total, lig_id, t_start)
        if tier == 1:
            state["seen"] += 1
            if state["seen"] >= 4:
                self.cancel()
        return result

    ScreeningRunner._emit_progress = emit_then_cancel
    try:
        report = runner.run(resume=False)
    finally:
        ScreeningRunner._emit_progress = original_emit

    accounting = runner.abort_accounting()
    aborted = set(accounting["aborted"])
    failed = set(accounting["failed"])
    assert aborted, "the cancel did not record any truncation"
    assert failed, "the engine failures were not recorded"
    assert not (aborted & failed), f"buckets overlap: {sorted(aborted & failed)}"
    assert report.is_conserved is False, _buckets(report, runner)

    total = report.n_funnel_filtered + report.n_dock_failed + len(aborted) + report.n_docked_success
    assert total == report.n_docking_pool, _buckets(report, runner)


def test_tiered_complete_run_without_incident_conserves(tmp_path: Path) -> None:
    """Reverse control: nothing failed and nothing was cancelled -> conserves."""
    engine = FaultEngine()
    runner = _runner(tmp_path, engine, 12)
    report = runner.run(resume=False)

    assert report.n_dock_failed == 0
    assert report.n_funnel_filtered == 0
    assert report.n_docked_success == 12
    assert report.is_conserved is True, _buckets(report, runner)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))
