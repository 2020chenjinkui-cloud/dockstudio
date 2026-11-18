"""Opt-in **real-machine** test for the gnina adapter (Linux + NVIDIA GPU only).

Why this file exists: until 2026-11-18 *every* gnina assertion in this repo was
either a fake ``#!/bin/sh`` binary (``test_gpu_adapters.py::test_gnina_mock_binary_e2e``)
or a CLI-assembly/parse check against hand-written fixtures
(``test_scoring_rigor.py::test_parse_gnina_*``).  The GPU code path had never
met a real gnina binary.  This module is that meeting point, and it is the only
place in the suite where a real binary, real CUDA and the adapter's own parser
are exercised together.

It is skipped unless ``DOCKSTUDIO_REAL_GPU=1``::

    DOCKSTUDIO_REAL_GPU=1 DOCKSTUDIO_GNINA=/path/to/gnina \\
        pytest server/tests/test_gnina_real_gpu.py -q

GN-04 (empty receptor pseudo-success) and GN-08 (in-place rescore) were fixed in
``fix/gnina-rescore``; those cases are plain assertions again.  What remains
``xfail(strict=False)`` is behaviour that is still open (GN-01 missing-file
exception type, GN-05 pose ordering) — see
``docs/qa/gnina-real-machine-2026-11-18.md`` for the measured evidence.
"""

from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent / "fixtures"

pytestmark = pytest.mark.skipif(
    os.environ.get("DOCKSTUDIO_REAL_GPU") != "1",
    reason="real-GPU test is opt-in: set DOCKSTUDIO_REAL_GPU=1 (needs Linux + NVIDIA GPU + gnina)",
)

BOX_SIZE = 22.5


def _gnina_binary() -> str | None:
    """Resolve gnina the way the adapter does, without touching CUDA."""
    env = os.environ.get("DOCKSTUDIO_GNINA")
    if env and Path(env).is_file():
        return env
    found = shutil.which("gnina")
    if found:
        return found
    # resources/bin lives at the *repo* root, next to server/ (config.py resolves
    # it the same way).  Looking only under server/ silently turned this whole
    # module into a no-op skip whenever DOCKSTUDIO_GNINA was not exported.
    test_file = Path(__file__).resolve()
    for root in (test_file.parents[2], test_file.parents[1]):
        for name in ("linux-x86_64", "linux-aarch64"):
            candidate = root / "resources" / "bin" / name / "gnina"
            if candidate.is_file():
                return str(candidate)
    return None


@pytest.fixture(scope="module")
def gnina_engine():
    if platform.system() != "Linux":
        pytest.skip("the gnina adapter is Linux-only (Windows users need WSL2)")
    if _gnina_binary() is None:
        pytest.skip("no gnina binary found (set DOCKSTUDIO_GNINA or resources/bin/linux-x86_64)")
    from dockstudio.engines.gnina import GninaEngine

    engine = GninaEngine()
    ok, reason = engine.is_available()
    if not ok:
        pytest.skip(f"gnina unavailable on this host: {reason}")
    return engine


@pytest.fixture(scope="module")
def inputs(tmp_path_factory) -> dict:
    """1HSG receptor + indinavir ligand, prepared through dockstudio's own code."""
    from dockstudio.prep.ligand import prepare_ligand_from_smiles
    from dockstudio.prep.receptor import prepare_receptor

    out_dir = tmp_path_factory.mktemp("gnina_real_inputs")
    receptor = prepare_receptor((FIXTURES / "1HSG.pdb").read_bytes(), "1HSG.pdb", out_dir)
    mk1 = next(c for c in receptor.cocrystal_ligands if c.resname == "MK1")
    ligand = prepare_ligand_from_smiles(
        (FIXTURES / "indinavir_smiles.txt").read_text(encoding="utf-8").strip(),
        "indinavir",
        out_dir,
    )
    return {
        "receptor": Path(receptor.pdbqt_path),
        "ligand": Path(ligand.pdbqt_path),
        "center": tuple(mk1.center),
    }


@pytest.fixture(scope="module")
def dock_run(gnina_engine, inputs, tmp_path_factory) -> dict:
    """One real docking run, shared by every assertion below (it costs ~30 s)."""
    from dockstudio.models import Box, DockParams

    work = tmp_path_factory.mktemp("gnina_real_dock")
    params = DockParams(exhaustiveness=2, num_modes=5, cpu=3, seed=42, scoring="vina")
    result = gnina_engine.dock(
        inputs["receptor"],
        inputs["ligand"],
        Box(center=inputs["center"], size=(BOX_SIZE, BOX_SIZE, BOX_SIZE)),
        params,
        work,
    )
    return {"result": result, "work": work}


# --------------------------------------------------------------------------
# The happy path: a real GPU docking run must produce real, sane numbers.
# --------------------------------------------------------------------------


def test_real_gpu_dock_returns_sane_poses(dock_run) -> None:
    result = dock_run["result"]
    work = dock_run["work"]

    assert result.engine == "gnina"
    assert result.engine_version not in ("", "unknown"), result.log[-2000:]
    assert result.engine_version.startswith("1."), result.engine_version
    assert (work / "out.pdbqt").is_file()

    # num_modes=5 asked for; a real 1HSG/indinavir run yields several modes
    assert len(result.poses) >= 2, result.log[-2000:]
    assert all(p.affinity < 0 for p in result.poses), [p.model_dump() for p in result.poses]
    # indinavir vs HIV-1 protease sits around -8..-13 kcal/mol
    assert -15.0 <= min(p.affinity for p in result.poses) <= -4.0, [
        p.model_dump() for p in result.poses
    ]
    # the banner + the mode table are the parts of stdout the adapter reports on
    assert "gnina v" in result.log
    assert "Commandline:" in result.log


def test_no_silent_cpu_fallback(dock_run) -> None:
    """A GPU host must not be downgraded to CPU behind the user's back.

    ``DOCKSTUDIO_GNINA_CPU`` is the only thing that may add ``--no_gpu``.
    """
    log = dock_run["result"].log
    cmd = next((ln for ln in log.splitlines() if ln.startswith("Commandline:")), "")
    assert cmd, log[-2000:]
    assert "--no_gpu" not in cmd, cmd
    # and the command the adapter actually built kept every mapped flag
    for flag in ("--receptor", "--ligand", "--center_x", "--size_x", "--exhaustiveness", "--out"):
        assert flag in cmd, cmd
    assert "--energy_range" not in cmd, "gnina has no --energy_range"


def test_cnn_scores_come_from_the_real_out_pdbqt(dock_run) -> None:
    """gnina's whole point: CNN pose score (0..1) and CNN affinity (pKd > 0).

    Real gnina 1.3.3 writes ``REMARK CNNscore`` / ``REMARK CNNaffinity`` into
    out.pdbqt, so the parser must find them there (the stdout table is a
    multi-line 5-column layout the fallback parser cannot be trusted with yet).
    """
    poses = dock_run["result"].poses
    assert poses[0].cnn_score is not None, "CNN score lost in parsing"
    assert poses[0].cnn_affinity is not None, "CNN affinity lost in parsing"
    assert 0.0 <= poses[0].cnn_score <= 1.0
    # pKd is positive and larger-is-better — NOT a kcal/mol energy (the mock
    # fixture in test_gpu_adapters.py uses -9.9 and gets the sign wrong).
    assert poses[0].cnn_affinity > 0.0


def test_adapter_parse_matches_the_output_file(dock_run) -> None:
    """Poses returned by dock() == an independent re-parse of the same file."""
    from dockstudio.engines.gnina import parse_gnina_out

    result = dock_run["result"]
    reparsed = parse_gnina_out(dock_run["work"] / "out.pdbqt", result.log)
    assert [p.index for p in reparsed] == [p.index for p in result.poses]
    assert [round(p.affinity, 4) for p in reparsed] == [round(p.affinity, 4) for p in result.poses]
    assert [p.cnn_score for p in reparsed] == [p.cnn_score for p in result.poses]


def test_exported_pose_keeps_coordinates(dock_run) -> None:
    """The exported SDF is what RMSD / pose QC consume: coordinates must survive."""
    from dockstudio.analysis.export import export_poses

    work = dock_run["work"]
    sdfs = export_poses(work / "out.pdbqt", work)
    assert len(sdfs) == len(dock_run["result"].poses)
    first_atom_pdbqt = next(
        ln for ln in (work / "out.pdbqt").read_text().splitlines() if ln.startswith("ATOM")
    )
    sdf_lines = sdfs[0].read_text().splitlines()
    counts_idx = next(i for i, ln in enumerate(sdf_lines) if "V2000" in ln)
    first_atom_sdf = sdf_lines[counts_idx + 1]
    assert [round(float(x), 3) for x in first_atom_pdbqt[30:54].split()] == [
        round(float(x), 3) for x in first_atom_sdf[0:30].split()
    ]


@pytest.mark.xfail(
    strict=False,
    reason="GN-05: gnina sorts poses by CNN score (--pose_sort_order=CNNscore), "
    "the adapter never re-sorts by affinity; every call site uses poses[0].affinity",
)
def test_poses_are_sorted_best_first_by_affinity(dock_run) -> None:
    affinities = [p.affinity for p in dock_run["result"].poses]
    assert affinities == sorted(affinities), affinities


# --------------------------------------------------------------------------
# Failure paths: the message a user sees must be honest and unambiguous.
# --------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=False,
    reason="GN-01: shutil.copy2's FileNotFoundError escapes dock() instead of DockingError",
)
def test_missing_ligand_raises_docking_error(gnina_engine, inputs, tmp_path: Path) -> None:
    from dockstudio.models import Box, DockingError, DockParams

    with pytest.raises(DockingError):
        gnina_engine.dock(
            inputs["receptor"],
            tmp_path / "does_not_exist.pdbqt",
            Box(center=inputs["center"], size=(BOX_SIZE, BOX_SIZE, BOX_SIZE)),
            DockParams(exhaustiveness=1, num_modes=1, cpu=1),
            tmp_path / "run",
        )


def test_empty_receptor_is_rejected(gnina_engine, inputs, tmp_path: Path) -> None:
    from dockstudio.models import Box, DockingError, DockParams

    empty = tmp_path / "empty_receptor.pdbqt"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(DockingError):
        gnina_engine.dock(
            empty,
            inputs["ligand"],
            Box(center=inputs["center"], size=(BOX_SIZE, BOX_SIZE, BOX_SIZE)),
            DockParams(exhaustiveness=1, num_modes=1, cpu=3),
            tmp_path / "run",
        )


def test_zero_size_box_raises_docking_error(gnina_engine, inputs, tmp_path: Path) -> None:
    """An invalid box must fail loudly — and the reason must be the real one."""
    from dockstudio.models import Box, DockingError, DockParams

    with pytest.raises(DockingError) as excinfo:
        gnina_engine.dock(
            inputs["receptor"],
            inputs["ligand"],
            Box(center=(0.0, 0.0, 0.0), size=(0.0, 0.0, 0.0)),
            DockParams(exhaustiveness=1, num_modes=1, cpu=3),
            tmp_path / "run",
        )
    msg = str(excinfo.value)
    assert "positive" in msg.lower(), msg
    # it must NOT be misreported as "box too large"
    assert "搜索盒设置超限" not in msg, msg


def test_in_place_rescore_accepts_an_exported_pose(gnina_engine, inputs, dock_run) -> None:
    """--score_only on the exported SDF — exactly what screen_runner Tier 2 does."""
    from dockstudio.analysis.export import export_poses
    from dockstudio.models import Box, DockParams

    work = dock_run["work"]
    pose1 = export_poses(work / "out.pdbqt", work)[0]
    scored = gnina_engine.dock(
        inputs["receptor"],
        pose1,
        Box(center=inputs["center"], size=(BOX_SIZE, BOX_SIZE, BOX_SIZE)),
        DockParams(exhaustiveness=1, num_modes=1, cpu=3, score_only=True),
        work / "score_only",
    )
    assert scored.poses, scored.log[-2000:]
    assert scored.poses[0].cnn_score is not None
    assert scored.poses[0].cnn_affinity is not None
    assert "GNINA_CNN_MISSING" not in scored.log
    assert "--score_only" in scored.log


# --------------------------------------------------------------------------
# DOCKSTUDIO_GNINA_CPU escape hatch (no GPU needed for these two checks).
# --------------------------------------------------------------------------


def test_cpu_opt_in_appends_no_gpu(gnina_engine, monkeypatch) -> None:
    from dockstudio.models import Box, DockParams

    monkeypatch.setenv("DOCKSTUDIO_GNINA_CPU", "1")
    cmd = gnina_engine._build_cmd(
        "gnina", Box(center=(0, 0, 0), size=(20, 20, 20)), DockParams(), "r", "l", "o"
    )
    assert "--no_gpu" in cmd


def test_cpu_opt_in_ungates_unavailable_gpu(gnina_engine, monkeypatch) -> None:
    """Without the opt-in a GPU-less host must say no; with it the engine must load."""
    monkeypatch.setattr(gnina_engine, "_gpu_ok", lambda: (False, "simulated: no NVIDIA GPU"))
    monkeypatch.delenv("DOCKSTUDIO_GNINA_CPU", raising=False)
    ok, reason = gnina_engine.is_available()
    assert not ok and reason, "a GPU-less host must not claim gnina is ready"

    monkeypatch.setenv("DOCKSTUDIO_GNINA_CPU", "1")
    ok, _ = gnina_engine.is_available()
    assert ok, "DOCKSTUDIO_GNINA_CPU=1 did not un-gate the engine"
