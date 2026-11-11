"""Adversarial stress-test suite for Milestone 3:
- Ligand path & filename collision stress (50 distinct ligands named ligand.pdbqt in separate folders)
- Partial failures & corrupted PDBQTs mixed into batch
- Edge cases: empty batch, single ligand, huge batch (500 ligands), argument-order polymorphism
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Any

import pytest

from dockstudio.engines.base import BaseEngine, EngineManifest
from dockstudio.engines.unidock import UniDockEngine
from dockstudio.models import (
    Box,
    DockingError,
    DockParams,
    DockResult,
    EngineUnavailableError,
    Pose,
    PreparedLigand,
)


_SAMPLE_PDBQT_TEMPLATE = """\
MODEL {model_idx}
REMARK VINA RESULT: {affinity:.2f} 0.000 0.000
ROOT
ATOM      1  C   UNL     1       {coord_x:.3f}   2.000   3.000  1.00  0.00    +0.000 C
ATOM      2  O   UNL     1       2.000   2.500   3.000  1.00  0.00    -0.300 O
ENDROOT
ENDMDL
"""


def make_sample_pdbqt(affinities: list[float], coord_offset: float = 1.0) -> str:
    blocks = []
    for idx, aff in enumerate(affinities, start=1):
        blocks.append(
            _SAMPLE_PDBQT_TEMPLATE.format(
                model_idx=idx, affinity=aff, coord_x=coord_offset + idx
            )
        )
    return "".join(blocks)


def create_fake_ligand(path: Path, affinities: list[float] | None = None, offset: float = 1.0) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    affs = affinities or [-8.0]
    path.write_text(make_sample_pdbqt(affs, coord_offset=offset), encoding="utf-8")
    return path


class MockSequentialEngine(BaseEngine):
    """Mock engine that implements dock() to test BaseEngine.dock_batch()."""

    manifest = EngineManifest(
        name="mock_seq",
        scoring_functions=["vina"],
        supports_batch_ligands=True,
    )

    def __init__(self, fail_ids: set[str] | None = None):
        self.fail_ids = fail_ids or set()
        self.docked: dict[str, Path] = {}

    def dock(
        self,
        receptor_pdbqt: Path,
        ligand_pdbqt: Path,
        box: Box,
        params: DockParams,
        work_dir: Path,
        flex_pdbqt: Path | None = None,
    ) -> DockResult:
        lig_id = ligand_pdbqt.stem
        if lig_id == "ligand":
            lig_id = ligand_pdbqt.parent.name

        if lig_id in self.fail_ids:
            raise DockingError(f"Simulated docking failure for {lig_id}")

        self.docked[lig_id] = work_dir
        # Parse affinity from input if available
        aff = -7.5
        if ligand_pdbqt.is_file():
            content = ligand_pdbqt.read_text(encoding="utf-8", errors="replace")
            for line in content.splitlines():
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    try:
                        aff = float(parts[3])
                        break
                    except (IndexError, ValueError):
                        pass

        out_file = work_dir / "out.pdbqt"
        out_file.write_text(make_sample_pdbqt([aff]), encoding="utf-8")
        return DockResult(
            engine="mock_seq",
            engine_version="1.0.0",
            params=params,
            box=box,
            poses=[Pose(index=1, affinity=aff)],
            log=f"docked {lig_id}",
            run_dir=str(work_dir),
        )


# =============================================================================
# Challenge 1: Ligand Path & Filename Collision Stress (50 ligands named ligand.pdbqt)
# =============================================================================

def test_unidock_collision_stress_50_ligands_all_named_ligand_pdbqt(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Stress-tests UniDockEngine.dock_batch with 50 distinct prepared ligands,
    where EVERY ligand file is named 'ligand.pdbqt' in a separate folder.
    Verifies:
    1. Zero overwrites in staging_ligands (50 unique staged files).
    2. Staged file content exactly matches respective source content.
    3. ligand_index.txt contains all 50 absolute staged paths without duplicates.
    4. dock_out results map accurately to each unique ligand_id with zero affinity swapping.
    """
    eng = UniDockEngine()
    rec_path = create_fake_ligand(tmp_path / "receptor.pdbqt")
    box = Box(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "work_collision_50"

    # Mock engine availability & version
    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    # 1. Create 50 distinct ligands: all named ligand.pdbqt in separate folders
    batch: list[tuple[str, Path]] = []
    expected_affinities: dict[str, float] = {}
    source_contents: dict[str, str] = {}

    for i in range(50):
        lig_id = f"mol_{i:03d}"
        aff = -5.0 - (i * 0.1)  # distinct affinities: -5.0, -5.1, ..., -9.9
        expected_affinities[lig_id] = round(aff, 2)
        lig_folder = tmp_path / "libraries" / "tcm_lib" / lig_id
        lig_file = lig_folder / "ligand.pdbqt"
        create_fake_ligand(lig_file, [expected_affinities[lig_id]], offset=float(i))
        source_contents[lig_id] = lig_file.read_text(encoding="utf-8")
        batch.append((lig_id, lig_file))

    # 2. Mock subprocess execution simulating Uni-Dock batch
    def fake_popen(cmd, *args, **kwargs):
        # Uni-Dock receives --ligand_index and --dir
        assert "--ligand_index" in cmd
        assert "--dir" in cmd
        idx_path = Path(cmd[cmd.index("--ligand_index") + 1])
        out_dir = Path(cmd[cmd.index("--dir") + 1])

        # Read the staged files from ligand_index
        lines = idx_path.read_text(encoding="utf-8").strip().splitlines()
        for line in lines:
            staged = Path(line.strip())
            stem = staged.stem  # e.g. mol_000
            # Read affinity from staged file
            staged_text = staged.read_text(encoding="utf-8")
            aff = -7.5
            for row in staged_text.splitlines():
                if "REMARK VINA RESULT:" in row:
                    aff = float(row.split()[3])
                    break
            # Uni-Dock outputs <staged_stem>_out.pdbqt
            (out_dir / f"{stem}_out.pdbqt").write_text(make_sample_pdbqt([aff]), encoding="utf-8")

        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return (f"Uni-Dock completed {len(lines)} ligands", "")

        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    # Verification 1: Staging directory contains exactly 50 files
    staging_dir = work_dir / "staging_ligands"
    assert staging_dir.is_dir()
    staged_files = list(staging_dir.glob("*.pdbqt"))
    assert len(staged_files) == 50, f"Expected 50 staged files, got {len(staged_files)}"

    # Verification 2: Zero overwrites & content fidelity
    for lig_id, _ in batch:
        staged_file = staging_dir / f"{lig_id}.pdbqt"
        assert staged_file.is_file(), f"Missing staged file for {lig_id}"
        # Staged file content must match source content exactly
        assert staged_file.read_text(encoding="utf-8") == source_contents[lig_id], (
            f"Content mismatch/overwrite in staged file for {lig_id}"
        )

    # Verification 3: ligand_index.txt integrity
    index_file = work_dir / "ligand_index.txt"
    assert index_file.is_file()
    index_lines = [l.strip() for l in index_file.read_text(encoding="utf-8").strip().splitlines() if l.strip()]
    assert len(index_lines) == 50, f"Expected 50 lines in ligand_index.txt, got {len(index_lines)}"
    # No duplicate lines
    assert len(set(index_lines)) == 50, "Duplicate paths found in ligand_index.txt"

    # Verification 4: Result mapping and zero affinity swapping
    assert len(results) == 50, f"Expected 50 results, got {len(results)}"
    for lig_id, exp_aff in expected_affinities.items():
        assert lig_id in results, f"Result missing for {lig_id}"
        actual_aff = results[lig_id].poses[0].affinity
        assert actual_aff == pytest.approx(exp_aff, abs=0.01), (
            f"Affinity swap detected for {lig_id}: expected {exp_aff}, got {actual_aff}"
        )


def test_base_engine_collision_stress_50_ligands_all_named_ligand_pdbqt(tmp_path: Path) -> None:
    """Stress-tests BaseEngine.dock_batch sequential fallback with 50 distinct prepared ligands
    where every file is named 'ligand.pdbqt' in a separate folder.
    Verifies that each ligand is docked inside an isolated sub-directory and zero overwrites occur.
    """
    eng = MockSequentialEngine()
    rec_path = create_fake_ligand(tmp_path / "receptor.pdbqt")
    box = Box(center=(0.0, 0.0, 0.0), size=(20.0, 20.0, 20.0))
    params = DockParams(exhaustiveness=8)
    work_dir = tmp_path / "base_collision_50"

    batch: list[tuple[str, Path]] = []
    expected_affinities: dict[str, float] = {}

    for i in range(50):
        lig_id = f"herb_cmpd_{i:03d}"
        aff = -6.0 - (i * 0.05)
        expected_affinities[lig_id] = round(aff, 2)
        lig_folder = tmp_path / "herbs" / lig_id
        lig_file = lig_folder / "ligand.pdbqt"
        create_fake_ligand(lig_file, [expected_affinities[lig_id]], offset=float(i))
        batch.append((lig_id, lig_file))

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    assert len(results) == 50
    # Verify 50 isolated directories
    for lig_id, exp_aff in expected_affinities.items():
        assert lig_id in results
        assert results[lig_id].poses[0].affinity == pytest.approx(exp_aff, abs=0.01)
        # Check sub_work_dir
        sub_dir = work_dir / f"dock_{lig_id}"
        assert sub_dir.is_dir()
        assert (sub_dir / "out.pdbqt").is_file()


def test_prepared_ligand_objects_in_unidock_batch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that dock_batch seamlessly handles list[PreparedLigand] objects,
    extracting id and pdbqt_path correctly for staging.
    """
    eng = UniDockEngine()
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "work_prep_objects"

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    prep_ligands: list[PreparedLigand] = []
    for i in range(10):
        p = create_fake_ligand(tmp_path / f"prep_{i}" / "ligand.pdbqt", [-7.0 - i])
        prep_ligands.append(
            PreparedLigand(
                id=f"tcm_{i}",
                source="file",
                pdbqt_path=str(p),
            )
        )

    def fake_popen(cmd, *args, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        for i in range(10):
            (out_dir / f"tcm_{i}_out.pdbqt").write_text(make_sample_pdbqt([-7.0 - i]), encoding="utf-8")
        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return ("Done", "")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(rec_path, prep_ligands, box, params, work_dir)
    assert len(results) == 10
    assert results["tcm_0"].poses[0].affinity == pytest.approx(-7.0)
    assert results["tcm_9"].poses[0].affinity == pytest.approx(-16.0)


# =============================================================================
# Challenge 2: Partial Failures & Corrupted Inputs Mixed into Batch
# =============================================================================

def test_unidock_partial_failures_mixed_corrupted_inputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Stress-tests UniDockEngine.dock_batch when valid ligands are mixed with:
    - Non-existent files
    - 0-byte empty files
    - Garbage binary / corrupted text
    - Ligands that crash / produce no output file
    - Ligands that produce unparseable corrupted poses
    Verifies that:
    1. Valid ligands dock successfully and are accurately parsed.
    2. Invalid / corrupted ligands are safely skipped without halting the batch.
    3. Even with non-zero process returncode, valid results are preserved and returned.
    """
    eng = UniDockEngine()
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "work_corrupted"

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    # Construct mixed batch
    # Valid ligands: v1, v2, v3
    v1 = create_fake_ligand(tmp_path / "v1.pdbqt", [-8.5])
    v2 = create_fake_ligand(tmp_path / "v2.pdbqt", [-7.2])
    v3 = create_fake_ligand(tmp_path / "v3.pdbqt", [-9.1])

    # Corrupted / invalid ligands:
    # 1. Non-existent file
    missing = tmp_path / "does_not_exist.pdbqt"
    # 2. 0-byte empty file
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("", encoding="utf-8")
    # 3. Binary garbage
    garbage = tmp_path / "garbage.pdbqt"
    garbage.write_bytes(b"\x00\xFF\xFE\x01\x02\x03BAD_HEADER\x00")
    # 4. Valid file that produces no output from unidock
    no_out = create_fake_ligand(tmp_path / "no_out.pdbqt", [-6.0])
    # 5. Valid file whose output file is unparseable garbage
    bad_out = create_fake_ligand(tmp_path / "bad_out.pdbqt", [-6.5])

    batch = [
        ("v1", v1),
        ("missing_file", missing),
        ("empty_file", empty),
        ("v2", v2),
        ("garbage_file", garbage),
        ("no_out_lig", no_out),
        ("bad_out_lig", bad_out),
        ("v3", v3),
    ]

    def fake_popen(cmd, *args, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        # Write valid outputs for v1, v2, v3
        (out_dir / "v1_out.pdbqt").write_text(make_sample_pdbqt([-8.5]), encoding="utf-8")
        (out_dir / "v2_out.pdbqt").write_text(make_sample_pdbqt([-7.2]), encoding="utf-8")
        (out_dir / "v3_out.pdbqt").write_text(make_sample_pdbqt([-9.1]), encoding="utf-8")

        # no_out_lig: write nothing
        # bad_out_lig: write corrupt unparseable text
        (out_dir / "bad_out_lig_out.pdbqt").write_text("CORRUPTED NOT A PDBQT AT ALL", encoding="utf-8")

        class FakeProcess:
            returncode = 1  # Non-zero exit code due to skipped/failed ligands!
            def communicate(self, timeout=None):
                return ("Uni-Dock: 3 succeeded, 3 failed", "Error parsing empty_file.pdbqt")

        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    # Missing file should not even be staged
    staging_dir = work_dir / "staging_ligands"
    assert not (staging_dir / "missing_file.pdbqt").exists()

    # Valid ligands MUST succeed
    assert "v1" in results
    assert results["v1"].poses[0].affinity == pytest.approx(-8.5)
    assert "v2" in results
    assert results["v2"].poses[0].affinity == pytest.approx(-7.2)
    assert "v3" in results
    assert results["v3"].poses[0].affinity == pytest.approx(-9.1)

    # Failed/corrupted ligands MUST NOT be in results
    assert "missing_file" not in results
    assert "no_out_lig" not in results
    assert "bad_out_lig" not in results

    assert len(results) == 3


def test_unidock_all_ligands_fail_raises_docking_error(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verifies that if ALL ligands in a batch fail and returncode != 0,
    UniDockEngine.dock_batch raises DockingError.
    """
    eng = UniDockEngine()
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "work_all_fail"

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    l1 = create_fake_ligand(tmp_path / "l1.pdbqt")
    batch = [("l1", l1)]

    def fake_popen(cmd, *args, **kwargs):
        class FakeProcess:
            returncode = 255
            def communicate(self, timeout=None):
                return ("", "CUDA out of memory error: failed to allocate grid")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    with pytest.raises(DockingError) as exc_info:
        eng.dock_batch(rec_path, batch, box, params, work_dir)

    assert "unidock batch exited with code 255" in str(exc_info.value)
    assert "CUDA" in str(exc_info.value) or "error" in str(exc_info.value).lower()


def test_base_engine_partial_failure_isolation(tmp_path: Path) -> None:
    """Verifies that BaseEngine.dock_batch isolates failures:
    failing individual ligands with simulated exceptions does not halt other ligands.
    """
    eng = MockSequentialEngine(fail_ids={"fail_1", "fail_2"})
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=8)
    work_dir = tmp_path / "base_partial_fail"

    good1 = create_fake_ligand(tmp_path / "good1" / "ligand.pdbqt", [-8.0])
    fail1 = create_fake_ligand(tmp_path / "fail_1" / "ligand.pdbqt", [-7.0])
    good2 = create_fake_ligand(tmp_path / "good2" / "ligand.pdbqt", [-9.0])
    fail2 = create_fake_ligand(tmp_path / "fail_2" / "ligand.pdbqt", [-6.0])
    missing = tmp_path / "nonexistent" / "ligand.pdbqt"

    batch = [
        ("good1", good1),
        ("fail_1", fail1),
        ("good2", good2),
        ("fail_2", fail2),
        ("missing", missing),
    ]

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    assert len(results) == 2
    assert "good1" in results
    assert "good2" in results
    assert results["good1"].poses[0].affinity == pytest.approx(-8.0)
    assert results["good2"].poses[0].affinity == pytest.approx(-9.0)
    assert "fail_1" not in results
    assert "fail_2" not in results
    assert "missing" not in results


# =============================================================================
# Challenge 3: Edge Case Inputs (Empty Batch, Single Ligand, Huge Batch)
# =============================================================================

def test_dock_batch_empty_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies both UniDockEngine and BaseEngine return empty dict immediately on empty batch."""
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)

    # 1. BaseEngine
    base_eng = MockSequentialEngine()
    res_base = base_eng.dock_batch(rec_path, [], box, params, tmp_path / "empty_base")
    assert res_base == {}

    # 2. UniDockEngine with binary mocked (simulating Linux)
    unidock_eng = UniDockEngine()
    monkeypatch.setattr(unidock_eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(unidock_eng, "_version", lambda binary: "1.1.0-mock")
    res_unidock = unidock_eng.dock_batch(rec_path, [], box, params, tmp_path / "empty_unidock")
    assert res_unidock == {}


def test_unidock_unsupported_platform_raises_unavailable(tmp_path: Path) -> None:
    """Verifies that calling UniDockEngine.dock_batch on Windows without WSL2
    cleanly raises EngineUnavailableError with helpful guidance.
    """
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    eng = UniDockEngine()

    if os.name == "nt":
        with pytest.raises(EngineUnavailableError) as exc_info:
            eng.dock_batch(rec_path, [("m1", rec_path)], box, params, tmp_path / "unavail")
        assert "Linux + NVIDIA GPU" in str(exc_info.value) or "WSL2" in str(exc_info.value)


def test_dock_batch_all_missing_files_input(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that if all items in batch point to non-existent files, empty dict is returned."""
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)

    batch = [("m1", tmp_path / "none1.pdbqt"), ("m2", tmp_path / "none2.pdbqt")]

    # BaseEngine
    base_eng = MockSequentialEngine()
    res_base = base_eng.dock_batch(rec_path, batch, box, params, tmp_path / "miss_base")
    assert res_base == {}

    # UniDockEngine
    unidock_eng = UniDockEngine()
    monkeypatch.setattr(unidock_eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(unidock_eng, "_version", lambda binary: "1.1.0-mock")
    res_unidock = unidock_eng.dock_batch(rec_path, batch, box, params, tmp_path / "miss_unidock")
    assert res_unidock == {}


def test_dock_batch_single_ligand(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies single-ligand batch works properly through dock_batch."""
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "single_run"

    eng = UniDockEngine()
    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    single = create_fake_ligand(tmp_path / "solo" / "ligand.pdbqt", [-10.5])
    batch = [("solo_cmpd", single)]

    def fake_popen(cmd, *args, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        (out_dir / "solo_cmpd_out.pdbqt").write_text(make_sample_pdbqt([-10.5]), encoding="utf-8")
        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return ("Done", "")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    res = eng.dock_batch(rec_path, batch, box, params, work_dir)
    assert len(res) == 1
    assert "solo_cmpd" in res
    assert res["solo_cmpd"].poses[0].affinity == pytest.approx(-10.5)

    index_lines = (work_dir / "ligand_index.txt").read_text(encoding="utf-8").strip().splitlines()
    assert len(index_lines) == 1
    assert "solo_cmpd.pdbqt" in index_lines[0]


def test_dock_batch_huge_batch_500_ligands_cli_length_and_timeout(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Stress-tests UniDockEngine.dock_batch with 500 ligands.
    Verifies:
    1. CLI command length is small (< 1000 chars, well below Windows 8191 limit)
       and constant O(1) complexity because --ligand_index is used.
    2. Dynamic timeout scales up to max(job_timeout_s, len(valid_items) * 10).
    3. ligand_index.txt contains exactly 500 lines.
    4. All 500 parsed results are returned.
    """
    eng = UniDockEngine()
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)
    work_dir = tmp_path / "huge_500_run"

    monkeypatch.setattr(eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(eng, "_version", lambda binary: "1.1.0-mock")

    batch: list[tuple[str, Path]] = []
    # Share a dummy file for disk efficiency, but with 500 distinct ligand IDs
    dummy_pdbqt = create_fake_ligand(tmp_path / "shared_dummy.pdbqt", [-8.0])
    for i in range(500):
        batch.append((f"huge_{i:04d}", dummy_pdbqt))

    captured_cmd: list[str] = []
    captured_timeout: list[float | None] = []

    def fake_popen(cmd, *args, **kwargs):
        captured_cmd.extend(cmd)
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        # Generate 500 mock outputs
        for i in range(500):
            (out_dir / f"huge_{i:04d}_out.pdbqt").write_text(make_sample_pdbqt([-8.0]), encoding="utf-8")

        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                captured_timeout.append(timeout)
                return ("Processed 500 ligands in batch", "")

        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    results = eng.dock_batch(rec_path, batch, box, params, work_dir)

    # 1. Verify CLI command length is well below Windows 8191 char limit
    cli_str = " ".join(captured_cmd)
    assert len(cli_str) < 1000, f"CLI command length too large ({len(cli_str)}): {cli_str}"
    assert "--ligand_index" in captured_cmd
    assert "--dir" in captured_cmd

    # 2. Verify dynamic timeout scaling: 500 * 10 = 5000 seconds
    assert len(captured_timeout) == 1
    assert captured_timeout[0] == 5000

    # 3. Verify ligand_index lines
    index_file = work_dir / "ligand_index.txt"
    index_lines = [l for l in index_file.read_text(encoding="utf-8").strip().splitlines() if l.strip()]
    assert len(index_lines) == 500

    # 4. Verify results
    assert len(results) == 500
    assert results["huge_0000"].poses[0].affinity == pytest.approx(-8.0)
    assert results["huge_0499"].poses[0].affinity == pytest.approx(-8.0)


def test_argument_order_polymorphism(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Verifies that dock_batch works whether called with:
    dock_batch(receptor, batch, ...) OR dock_batch(batch, receptor, ...)
    on both UniDockEngine and BaseEngine.
    """
    rec_path = create_fake_ligand(tmp_path / "rec.pdbqt")
    lig_path = create_fake_ligand(tmp_path / "lig.pdbqt", [-7.5])
    box = Box(center=(0, 0, 0), size=(20, 20, 20))
    params = DockParams(exhaustiveness=4)

    # 1. BaseEngine
    base_eng = MockSequentialEngine()
    # Normal order: (rec, batch)
    r1 = base_eng.dock_batch(rec_path, [("m1", lig_path)], box, params, tmp_path / "b1")
    assert "m1" in r1
    # Inverted order: (batch, rec)
    r2 = base_eng.dock_batch([("m1", lig_path)], rec_path, box, params, tmp_path / "b2")
    assert "m1" in r2

    # 2. UniDockEngine
    unidock_eng = UniDockEngine()
    monkeypatch.setattr(unidock_eng, "_binary", lambda: Path("/bin/unidock"))
    monkeypatch.setattr(unidock_eng, "_version", lambda binary: "1.1.0-mock")

    def fake_popen(cmd, *args, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        (out_dir / "m1_out.pdbqt").write_text(make_sample_pdbqt([-7.5]), encoding="utf-8")
        class FakeProcess:
            returncode = 0
            def communicate(self, timeout=None):
                return ("Done", "")
        return FakeProcess()

    monkeypatch.setattr("subprocess.Popen", fake_popen)

    # Normal order
    u1 = unidock_eng.dock_batch(rec_path, [("m1", lig_path)], box, params, tmp_path / "u1")
    assert "m1" in u1
    # Inverted order
    u2 = unidock_eng.dock_batch([("m1", lig_path)], rec_path, box, params, tmp_path / "u2")
    assert "m1" in u2
