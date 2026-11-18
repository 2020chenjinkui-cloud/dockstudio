"""Tests for the Uni-Dock batch receptor guard (T9-F2).

``VinaEngine.dock`` validates the receptor, but ``UniDockEngine.dock_batch`` is a
separate native GPU entry point that does not go through it.  An atom-less
receptor makes Uni-Dock "succeed" with affinity 0.0 for every ligand, which
downstream cannot tell apart from a real hit (GN-04).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dockstudio.engines.unidock import UniDockEngine
from dockstudio.models import Box, DockingError, DockParams

GOOD_RECEPTOR = (
    "ATOM      1  N   ILE A  10      13.000  22.000   5.000  1.00  0.00    -0.300 N\n"
    "ATOM      2  CA  ILE A  10      14.000  22.500   5.200  1.00  0.00     0.100 C\n"
)


def _batch_dirs(tmp_path: Path):
    lig = tmp_path / "ligand.pdbqt"
    lig.write_text("ROOT\nENDROOT\n", encoding="utf-8")
    work = tmp_path / "run"
    work.mkdir()
    return [(("LIG_001"), lig)], work


def test_unidock_batch_rejects_empty_receptor(tmp_path: Path) -> None:
    """A 0-byte receptor must fail loudly instead of scoring every ligand 0.0."""
    batch, work = _batch_dirs(tmp_path)
    empty = tmp_path / "empty.pdbqt"
    empty.write_text("", encoding="utf-8")

    eng = UniDockEngine()
    with pytest.raises(DockingError, match="原子"):
        eng.dock_batch(
            empty,
            batch,
            Box(center=(0, 0, 0), size=(20, 20, 20)),
            DockParams(exhaustiveness=2, num_modes=1, seed=42),
            work,
        )


def test_unidock_batch_rejects_missing_receptor(tmp_path: Path) -> None:
    """A non-existent receptor path is a clear DockingError, not a crash."""
    batch, work = _batch_dirs(tmp_path)
    eng = UniDockEngine()
    with pytest.raises(DockingError):
        eng.dock_batch(
            tmp_path / "nope.pdbqt",
            batch,
            Box(center=(0, 0, 0), size=(20, 20, 20)),
            DockParams(exhaustiveness=2, num_modes=1, seed=42),
            work,
        )


def test_unidock_batch_accepts_receptor_with_atoms(tmp_path: Path) -> None:
    """The guard must not reject a valid receptor (fails later on the GPU, not here)."""
    batch, work = _batch_dirs(tmp_path)
    rec = tmp_path / "rec.pdbqt"
    rec.write_text(GOOD_RECEPTOR, encoding="utf-8")

    eng = UniDockEngine()
    # The binary is absent on this host, so the call stops with an engine error -
    # the point is that it got *past* the receptor guard (no "原子" message).
    with pytest.raises(Exception) as excinfo:
        eng.dock_batch(
            rec,
            batch,
            Box(center=(0, 0, 0), size=(20, 20, 20)),
            DockParams(exhaustiveness=2, num_modes=1, seed=42),
            work,
        )
    assert "原子" not in str(excinfo.value)
