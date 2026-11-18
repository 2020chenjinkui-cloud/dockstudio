"""GPU-free regressions for the GNINA in-place rescore path (GN-08/GN-03/GN-04).

Everything here runs on a CPU-only machine: the "gnina" binary is a shell stub
that reproduces the *input validation* of the real gnina 1.3.3, so the failure
mode these tests reproduce is the one measured on the P4:

    Parse error on line 2 in file "ligand.pdbqt": Unknown or inappropriate tag

Before the fix, ``VinaEngine.dock()`` renamed every ligand to ``ligand.pdbqt``,
gnina's PDBQT reader parsed the exported ``pose_1.sdf`` as PDBQT, and
``screen_runner`` swallowed the resulting ``DockingError`` with
``except Exception: pass`` — CNN columns stayed empty with no trace anywhere.

Each test in this module is written to fail on the pre-fix code:
 * the SDF-rescore test fails with the real ``Parse error``,
 * the multi-model test fails on ``ENDMDL``,
 * the degradation-marker tests fail because ``gnina_cnn_rescore`` did not exist,
 * the empty-receptor test fails because a pseudo pose was returned.
"""

from __future__ import annotations

import platform
import sys
import textwrap
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockstudio.core import screen_runner as sr
from dockstudio.engines.base import (
    DEGRADED_PREFIX,
    GNINA_CNN_MISSING,
    GNINA_RESCORE_UNAVAILABLE,
    marker,
)
from dockstudio.engines.gnina import GninaEngine
from dockstudio.models import (
    Box,
    DockingError,
    DockParams,
    DockResult,
    Pose,
    PreparedLigand,
    ScreeningCandidateRow,
    ScreeningConfig,
)

BOX = Box(center=(13.073, 22.467, 5.557), size=(22.5, 22.5, 22.5))

RECEPTOR_PDBQT = textwrap.dedent(
    """\
    ATOM      1  N   ILE A  10      13.000  22.000   5.000  1.00  0.00    -0.300 N
    ATOM      2  CA  ILE A  10      14.000  22.500   5.200  1.00  0.00     0.100 C
    END
    """
)

#: A genuine (if tiny) SDF, exactly the shape ``export_poses`` writes for a pose.
SDF_POSE = textwrap.dedent(
    """\

         RDKit          3D

      3  2  0  0  0  0  0  0  0  0999 V2000
       14.1730   24.0660    4.4220 C   0  0  0  0  0  0  0  0  0  0  0  0
       12.7950   23.5960    4.9380 C   0  0  0  0  0  0  0  0  0  0  0  0
       11.6870   24.5540    4.6420 N   0  0  0  0  0  0  0  0  0  0  0  0
      1  2  1  0
      2  3  1  0
    M  END
    """
)

#: ``out.pdbqt`` from a real run: two MODEL blocks, each closed by ENDMDL.
MULTI_MODEL_PDBQT = textwrap.dedent(
    """\
    MODEL 1
    REMARK minimizedAffinity -11.194829
    REMARK CNNscore 0.860957563
    REMARK CNNaffinity 8.77859974
    ROOT
    ATOM      1  C   UNL     1      14.173  24.066   4.422  1.00  0.00    +0.129 C
    ENDROOT
    ENDMDL
    MODEL 2
    REMARK minimizedAffinity -10.4713345
    REMARK CNNscore 0.737501442
    REMARK CNNaffinity 9.17348194
    ROOT
    ATOM      1  C   UNL     1      15.100  23.100   5.100  1.00  0.00    +0.129 C
    ENDROOT
    ENDMDL
    """
)


_FAKE_GNINA_PY = """#!__INTERPRETER__
# Cross-platform stand-in for gnina 1.3.3 (T9-F1): a POSIX ``#!/bin/sh`` stub
# cannot run on the windows-latest leg of ci.yml's server-test-full matrix.
import sys
from pathlib import Path

CNN = __CNN__
args = sys.argv[1:]
if "--version" in args:
    print("gnina 1.3.3")
    raise SystemExit(0)

Path("argv.txt").write_text("\\n".join(args) + "\\n", encoding="utf-8")
prev = ""
out = lig = ""
for a in args:
    if prev == "--out":
        out = a
    if prev == "--ligand":
        lig = a
    prev = a

if lig.endswith(".pdbqt"):
    text = Path(lig).read_text(encoding="utf-8", errors="replace")
    if "V2000" in text:
        sys.stderr.write(f'Parse error on line 2 in file "{lig}": Unknown or inappropriate tag\\n')
        raise SystemExit(1)
    if any(line.startswith("ENDMDL") for line in text.splitlines()):
        sys.stderr.write(f'Parse error on line 92 in file "{lig}": Unknown or inappropriate tag\\n')
        raise SystemExit(1)

lines = ["MODEL 1", "REMARK minimizedAffinity -11.194829"]
if CNN:
    lines += ["REMARK CNNscore 0.860957563", "REMARK CNNaffinity 8.77859974"]
lines += ["ROOT", "ENDROOT", "ENDMDL"]
Path(out).write_text("\\n".join(lines) + "\\n", encoding="utf-8")
print("gnina v1.3.3 (stub)")
"""


# T12-F1: Windows ``CreateProcess`` has no shell and ignores the shebang, so a
# ``#!<python>`` script cannot be launched as a program - measured on Linux by
# stripping the shebang (the kernel answers ENOEXEC, exactly like Windows), which
# makes the five stub-executing tests below fail.  They are skipped there instead
# of pretending the stub is portable.  The staging / guard / parsing logic that
# does *not* spawn a process keeps running on every platform.
_NEEDS_EXECUTABLE_STUB = pytest.mark.skipif(
    sys.platform == "win32",
    reason=(
        "T12-F1: the gnina stub is a #!<python> script; Windows CreateProcess cannot "
        "execute it (no shell, shebang ignored). Covered on Linux/macOS instead."
    ),
)


def _fake_gnina(tmp_path: Path, *, cnn: bool = True, name: str = "fake_gnina") -> Path:
    """Cross-platform stand-in for the gnina binary that validates its input.

    Implemented as a Python script with an ``env python`` shebang rather than a
    ``/bin/sh`` script: Windows has no ``sh``, so the old POSIX stub made the
    whole ``server-test-full`` matrix red on windows-latest (finding T9-F1).
    The behaviour under test - argv handling, PDBQT-vs-SDF rejection, CNN
    REMARKs - is unchanged.

    ``cnn=False`` makes it emit an ``out.pdbqt`` without any CNN REMARK, which is
    how a degraded run looks (GN-03).
    """
    script = tmp_path / name
    script.write_text(
        _FAKE_GNINA_PY.replace("__CNN__", repr(cnn)).replace("__INTERPRETER__", sys.executable),
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


@pytest.fixture()
def staged_engine(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """GninaEngine wired to the stub binary on a pretend-Linux host."""
    binary = _fake_gnina(tmp_path)
    monkeypatch.setenv("DOCKSTUDIO_GNINA", str(binary))
    monkeypatch.delenv("DOCKSTUDIO_GNINA_CPU", raising=False)
    monkeypatch.setattr(platform, "system", lambda: "Linux")
    engine = GninaEngine()
    monkeypatch.setattr(engine, "_gpu_ok", lambda: (True, ""))
    return engine


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(RECEPTOR_PDBQT, encoding="utf-8")
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text("ROOT\nENDROOT\n", encoding="utf-8")
    return receptor, ligand


# ---------------------------------------------------------------------------
# GN-08: format handling. These are the tests that were RED before the fix.
# ---------------------------------------------------------------------------


def test_sdf_pose_keeps_its_extension_when_staged(tmp_path: Path) -> None:
    """An exported pose must not be renamed into ``ligand.pdbqt``."""
    src = tmp_path / "pose_1.sdf"
    src.write_text(SDF_POSE, encoding="utf-8")
    work = tmp_path / "run"
    work.mkdir()
    staged = GninaEngine()._stage_ligand(src, work)
    assert staged.name == "ligand.sdf"
    assert staged.read_text(encoding="utf-8") == SDF_POSE


@_NEEDS_EXECUTABLE_STUB
def test_in_place_rescore_accepts_an_exported_sdf_pose(
    staged_engine: GninaEngine, tmp_path: Path
) -> None:
    """The Tier-2 call shape: ``--score_only`` on ``pose_1.sdf`` must not fail."""
    receptor, _ = _inputs(tmp_path)
    pose = tmp_path / "pose_1.sdf"
    pose.write_text(SDF_POSE, encoding="utf-8")

    result = staged_engine.dock(
        receptor,
        pose,
        BOX,
        DockParams(exhaustiveness=1, num_modes=1, score_only=True),
        tmp_path / "run",
    )

    assert result.poses, result.log[-2000:]
    assert result.poses[0].cnn_score == pytest.approx(0.860957563)
    assert result.poses[0].cnn_affinity == pytest.approx(8.77859974)
    # the subprocess really received the in-place flag *and* the SDF file name
    argv = (tmp_path / "run" / "argv.txt").read_text(encoding="utf-8").split()
    assert "--score_only" in argv, argv
    assert "ligand.sdf" in argv, argv
    assert "ligand.pdbqt" not in argv, argv


@_NEEDS_EXECUTABLE_STUB
def test_multimodel_out_pdbqt_is_reduced_to_one_model(
    staged_engine: GninaEngine, tmp_path: Path
) -> None:
    """A multi-model ``out.pdbqt`` must be rescorable (MODEL 1, no ENDMDL)."""
    receptor, _ = _inputs(tmp_path)
    rescore_input = tmp_path / "out.pdbqt"
    rescore_input.write_text(MULTI_MODEL_PDBQT, encoding="utf-8")

    result = staged_engine.dock(
        receptor,
        rescore_input,
        BOX,
        DockParams(exhaustiveness=1, num_modes=1, minimize=True),
        tmp_path / "run",
    )

    staged = (tmp_path / "run" / "ligand.pdbqt").read_text(encoding="utf-8")
    assert staged.count("MODEL") == 1, staged
    assert "ENDMDL" not in staged, staged
    # MODEL 1 is the documented choice: its minimised affinity is what stays
    assert "REMARK minimizedAffinity -11.194829" in staged
    assert "REMARK minimizedAffinity -10.4713345" not in staged
    assert result.poses[0].affinity == pytest.approx(-11.194829)


def test_plain_single_pose_pdbqt_is_left_untouched(tmp_path: Path) -> None:
    """The normal meeko ligand must not be rewritten by the normaliser."""
    plain = "ROOT\nATOM      1  C   UNL     1       1.000   2.000   3.000\nENDROOT\nTORSDOF 0\n"
    src = tmp_path / "ligand.pdbqt"
    src.write_text(plain, encoding="utf-8")
    work = tmp_path / "run"
    work.mkdir()
    staged = GninaEngine()._stage_ligand(src, work)
    assert staged.read_text(encoding="utf-8") == plain


# ---------------------------------------------------------------------------
# GN-04: an atom-less receptor must be rejected, not "docked"
# ---------------------------------------------------------------------------


def test_staging_never_rewrites_the_callers_file(tmp_path: Path) -> None:
    """A retry that points at our own staging name must not truncate the source."""
    work = tmp_path / "run"
    work.mkdir()
    source = work / "ligand.pdbqt"
    source.write_text(MULTI_MODEL_PDBQT, encoding="utf-8")

    staged = GninaEngine()._stage_ligand(source, work)

    assert staged.name == "ligand_staged.pdbqt"
    assert source.read_text(encoding="utf-8") == MULTI_MODEL_PDBQT  # untouched
    assert staged.read_text(encoding="utf-8").count("MODEL") == 1


@_NEEDS_EXECUTABLE_STUB
def test_empty_receptor_is_rejected_before_launching(
    staged_engine: GninaEngine, tmp_path: Path
) -> None:
    """0-byte receptor used to yield affinity == 0.0 with plausible CNN values."""
    _, ligand = _inputs(tmp_path)
    empty = tmp_path / "empty_receptor.pdbqt"
    empty.write_text("", encoding="utf-8")

    with pytest.raises(DockingError) as excinfo:
        staged_engine.dock(
            empty,
            ligand,
            BOX,
            DockParams(exhaustiveness=1, num_modes=1),
            tmp_path / "run",
        )

    message = str(excinfo.value)
    assert "不含任何原子记录" in message, message
    assert "空受体无法用于对接" in message, message
    # the stub never ran: no out.pdbqt and no staged receptor copy
    assert not (tmp_path / "run" / "out.pdbqt").exists()


@_NEEDS_EXECUTABLE_STUB
def test_receptor_with_atoms_is_accepted(staged_engine: GninaEngine, tmp_path: Path) -> None:
    """The guard must not over-reject a usable receptor."""
    receptor, ligand = _inputs(tmp_path)
    result = staged_engine.dock(
        receptor, ligand, BOX, DockParams(exhaustiveness=1, num_modes=1), tmp_path / "run"
    )
    assert result.poses


# ---------------------------------------------------------------------------
# GN-03: missing CNN values must be announced, not silently None
# ---------------------------------------------------------------------------


@_NEEDS_EXECUTABLE_STUB
def test_cnn_missing_result_carries_a_degradation_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    binary = _fake_gnina(tmp_path, cnn=False, name="fake_gnina_no_cnn")
    monkeypatch.setenv("DOCKSTUDIO_GNINA", str(binary))
    monkeypatch.setattr(platform, "system", lambda: "Linux")
    engine = GninaEngine()
    monkeypatch.setattr(engine, "_gpu_ok", lambda: (True, ""))
    receptor, ligand = _inputs(tmp_path)

    result = engine.dock(
        receptor, ligand, BOX, DockParams(exhaustiveness=1, num_modes=1), tmp_path / "run"
    )

    assert result.poses[0].cnn_score is None
    assert f"{DEGRADED_PREFIX} {GNINA_CNN_MISSING}" in result.log, result.log[-1500:]


# ---------------------------------------------------------------------------
# The Tier-3 rescore helper: markers instead of `except Exception: pass`
# ---------------------------------------------------------------------------


class _StubGnina:
    """Minimal DockEngine stand-in that returns/raises on demand."""

    def __init__(self, *, error: Exception | None = None, poses: list[Pose] | None = None):
        self.error = error
        self.poses = poses if poses is not None else []
        self.calls = 0

    def dock(self, *_args, **_kwargs) -> DockResult:
        self.calls += 1
        if self.error is not None:
            raise self.error
        return DockResult(
            engine="gnina",
            engine_version="1.3.3",
            params=DockParams(),
            box=BOX,
            poses=self.poses,
            log="",
            run_dir="/tmp/gnina_stub",
        )


def _pose(**kwargs) -> Pose:
    base = {"index": 1, "affinity": -11.19, "cnn_score": 0.86, "cnn_affinity": 8.78}
    base.update(kwargs)
    return Pose(**base)


def test_rescore_failure_returns_an_unavailable_marker(tmp_path: Path) -> None:
    """The exact production failure: a parse error that used to be swallowed."""
    pose_file = tmp_path / "pose_1.sdf"
    pose_file.write_text(SDF_POSE, encoding="utf-8")
    engine = _StubGnina(
        error=DockingError('vina exited with code 1: Parse error on line 2 in file "ligand.pdbqt"')
    )

    cnn_score, cnn_affinity, degraded = sr.gnina_cnn_rescore(
        engine, tmp_path / "receptor.pdbqt", "lig1", pose_file, BOX, tmp_path / "run"
    )

    assert (cnn_score, cnn_affinity) == (None, None)
    assert degraded is not None
    assert degraded.startswith(f"{GNINA_RESCORE_UNAVAILABLE}:"), degraded
    assert "Parse error" in degraded
    assert engine.calls == 1


def test_rescore_missing_pose_file_returns_an_unavailable_marker(tmp_path: Path) -> None:
    engine = _StubGnina()
    _, _, degraded = sr.gnina_cnn_rescore(
        engine, tmp_path / "receptor.pdbqt", "lig1", None, BOX, tmp_path / "run"
    )
    assert degraded is not None
    assert degraded.startswith(f"{GNINA_RESCORE_UNAVAILABLE}:")
    assert "没有可用的对接构象文件" in degraded
    assert engine.calls == 0


def test_rescore_without_cnn_values_returns_a_cnn_missing_marker(tmp_path: Path) -> None:
    pose_file = tmp_path / "pose_1.sdf"
    pose_file.write_text(SDF_POSE, encoding="utf-8")
    engine = _StubGnina(poses=[_pose(cnn_score=None, cnn_affinity=None)])

    cnn_score, cnn_affinity, degraded = sr.gnina_cnn_rescore(
        engine, tmp_path / "receptor.pdbqt", "lig1", pose_file, BOX, tmp_path / "run"
    )

    assert (cnn_score, cnn_affinity) == (None, None)
    assert degraded is not None
    assert degraded.startswith(f"{GNINA_CNN_MISSING}:"), degraded


def test_rescore_success_returns_values_and_no_marker(tmp_path: Path) -> None:
    pose_file = tmp_path / "pose_1.sdf"
    pose_file.write_text(SDF_POSE, encoding="utf-8")
    engine = _StubGnina(poses=[_pose()])

    cnn_score, cnn_affinity, degraded = sr.gnina_cnn_rescore(
        engine, tmp_path / "receptor.pdbqt", "lig1", pose_file, BOX, tmp_path / "run"
    )

    assert (cnn_score, cnn_affinity) == (0.86, 8.78)
    assert degraded is None


# ---------------------------------------------------------------------------
# Landing points: candidate row alerts + screening_results.csv "Alerts" column
# ---------------------------------------------------------------------------


def _tier2_item(tmp_path: Path, ligand_id: str = "lig1") -> tuple[dict, PreparedLigand]:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    pose_path = tmp_path / f"pose_{ligand_id}.sdf"
    pose_path.write_text(Chem.MolToMolBlock(mol), encoding="utf-8")
    ligand = PreparedLigand(
        id=ligand_id,
        source="smiles",
        pdbqt_path=str(tmp_path / f"{ligand_id}.pdbqt"),
        name=ligand_id,
        smiles="CCO",
        min_energy_kcal=None,
    )
    (tmp_path / f"{ligand_id}.pdbqt").write_text("ROOT\nENDROOT\n", encoding="utf-8")
    item = {
        "ligand": ligand,
        "affinity_mean": -7.5,
        "affinity_std": 0.2,
        "affinity_seed0": -7.4,
        "pose_rmsd_seed": 0.4,
        "run_dir": str(tmp_path / "tier2_run"),
        "best_pose_path": str(pose_path),
        "sile": 1.0,
        "le": 0.3,
        "lle": 2.0,
    }
    return item, ligand


def test_tier3_candidate_carries_the_rescore_degradation_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end through ``_run_tier3``: the swallowed exception must surface."""
    item, ligand = _tier2_item(tmp_path)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(RECEPTOR_PDBQT, encoding="utf-8")
    stub_gnina = _StubGnina(error=DockingError("vina exited with code 1: Parse error on line 2"))
    monkeypatch.setattr(sr, "get_engine", lambda name: stub_gnina)
    monkeypatch.setattr(stub_gnina, "is_available", lambda: (True, None), raising=False)

    runner = sr.ScreeningRunner(
        config=ScreeningConfig(enable_gnina_rescore=True, box=BOX),
        receptor_pdbqt=receptor,
        ligands=[ligand],
        out_dir=tmp_path / "out",
        engine=stub_gnina,
    )
    candidates = runner._run_tier3([item], {})

    assert len(candidates) == 1
    alerts = candidates[0].alerts
    assert any(a.startswith(f"{GNINA_RESCORE_UNAVAILABLE}:") for a in alerts), alerts
    assert candidates[0].cnn_score is None

    csv_path = sr.export_screening_results_csv(candidates, tmp_path / "screening_results.csv")
    csv_text = csv_path.read_text(encoding="utf-8")
    assert "Alerts" in csv_text.splitlines()[0]
    assert GNINA_RESCORE_UNAVAILABLE in csv_text, csv_text.splitlines()[0]


def test_tier3_candidate_carries_cnn_missing_marker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    item, ligand = _tier2_item(tmp_path)
    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(RECEPTOR_PDBQT, encoding="utf-8")
    stub_gnina = _StubGnina(poses=[_pose(cnn_score=None, cnn_affinity=None)])
    monkeypatch.setattr(sr, "get_engine", lambda name: stub_gnina)
    monkeypatch.setattr(stub_gnina, "is_available", lambda: (True, None), raising=False)

    runner = sr.ScreeningRunner(
        config=ScreeningConfig(enable_gnina_rescore=True, box=BOX),
        receptor_pdbqt=receptor,
        ligands=[ligand],
        out_dir=tmp_path / "out",
        engine=stub_gnina,
    )
    candidates = runner._run_tier3([item], {})

    assert any(a.startswith(f"{GNINA_CNN_MISSING}:") for a in candidates[0].alerts), candidates[
        0
    ].alerts


def test_marker_helper_is_csv_safe() -> None:
    text = marker(GNINA_RESCORE_UNAVAILABLE, "line1\nline2 , with, commas; and semicolons")
    assert "\n" not in text
    assert text.startswith(f"{GNINA_RESCORE_UNAVAILABLE}:line1 line2")
    assert len(marker(GNINA_RESCORE_UNAVAILABLE, "x" * 1000)) < 300


def test_candidate_row_keeps_alert_roundtrip() -> None:
    """Alerts survive the model round-trip used by the JSONL checkpoint."""
    row = ScreeningCandidateRow(ligand_id="lig1", alerts=[marker(GNINA_CNN_MISSING, "no cnn")])
    again = ScreeningCandidateRow(**row.model_dump())
    assert again.alerts == row.alerts
