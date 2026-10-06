"""Unit tests for the vina PDBQT result parser."""

from __future__ import annotations

from pathlib import Path

import pytest

from dockstudio.engines.vina import parse_vina_out
from dockstudio.models import DockingError

SAMPLE = """\
MODEL 1
REMARK VINA RESULT:    -9.8      0.000      0.000
ROOT
ATOM      1  C   UNL     1       1.0 2.0 3.0
ENDROOT
ENDMDL
MODEL 2
REMARK VINA RESULT:    -9.1      2.031      4.502
ROOT
ATOM      1  C   UNL     1       1.0 2.0 3.0
ENDROOT
ENDMDL
MODEL 3
REMARK VINA RESULT:    -8.7      1.500      9.001
ROOT
ATOM      1  C   UNL     1       1.0 2.0 3.0
ENDROOT
ENDMDL
"""


def test_parse_three_models(tmp_path: Path) -> None:
    out = tmp_path / "out.pdbqt"
    out.write_text(SAMPLE, encoding="utf-8")
    poses = parse_vina_out(out)
    assert len(poses) == 3
    assert poses[0].index == 1 and poses[0].affinity == pytest.approx(-9.8)
    assert poses[0].rmsd_lb == pytest.approx(0.0)
    assert poses[2].index == 3 and poses[2].affinity == pytest.approx(-8.7)
    assert poses[2].rmsd_ub == pytest.approx(9.001)


def test_no_poses_raises(tmp_path: Path) -> None:
    out = tmp_path / "out.pdbqt"
    out.write_text("garbage\n", encoding="utf-8")
    with pytest.raises(DockingError):
        parse_vina_out(out)
