"""Unit tests for the HTML report generator."""

from __future__ import annotations

from datetime import datetime
from html.parser import HTMLParser
from pathlib import Path

from dockstudio.analysis.report import generate_report
from dockstudio.models import (
    Box,
    DockJobRequest,
    DockParams,
    DockResult,
    Interaction,
    Pose,
    PoseAnalysis,
    Project,
)


class _Checker(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags: set[str] = set()
        self.text = ""

    def handle_starttag(self, tag, attrs):
        self.tags.add(tag)

    def handle_data(self, data):
        self.text += data


def _sample_result() -> DockResult:
    return DockResult(
        engine="vina",
        engine_version="1.2.7",
        params=DockParams(exhaustiveness=4, num_modes=5, seed=42),
        box=Box(center=(13.0, 22.4, 5.5), size=(22.5, 22.5, 22.5)),
        poses=[
            Pose(index=1, affinity=-9.85, rmsd_lb=0.0, rmsd_ub=0.0),
            Pose(index=2, affinity=-9.12, rmsd_lb=2.1, rmsd_ub=8.4),
        ],
        run_dir="runs/x",
    )


def test_report_contains_key_fields(tmp_path: Path) -> None:
    receptor_pdb = tmp_path / "rec.pdb"
    receptor_pdb.write_text(
        "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    pose_sdf = tmp_path / "pose_1.sdf"
    pose_sdf.write_text("test\n  DockStudio\n\n  0  0  0\nM  END\n$$$$\n", encoding="utf-8")

    analysis = PoseAnalysis(
        pose_index=1,
        source="builtin",
        interactions=[
            Interaction(type="hbond", residue="ASP25:A", distance=2.9, detail="蛋白供体"),
            Interaction(type="hydrophobic", residue="ILE50:B", distance=3.6),
        ],
    )
    project = Project(id="p1", name="HIV 蛋白酶", path=str(tmp_path), created_at=datetime.now())
    request = DockJobRequest(
        receptor_id="rec123",
        ligand_id="lig456",
        box=Box(center=(13.0, 22.4, 5.5), size=(22.5, 22.5, 22.5)),
    )
    page = generate_report(
        project=project,
        request=request,
        result=_sample_result(),
        pose_analyses=[analysis],
        receptor_pdb_path=receptor_pdb,
        pose_sdf_paths=[pose_sdf],
        receptor_name="1HSG",
        ligand_name="indinavir",
    )
    assert page.startswith("<!DOCTYPE html>")
    for needle in (
        "HIV 蛋白酶",
        "vina 1.2.7",
        "-9.85",
        "ASP25:A",
        "氢键",
        "3Dmol",
        "免责声明",
        "indinavir",
        "rec123",
        "lig456",
    ):
        assert needle in page, f"missing {needle!r} in report"

    # the HTML must be parseable
    checker = _Checker()
    checker.feed(page)
    assert {"html", "table", "script", "div"} <= checker.tags
    assert "分子对接报告" in checker.text

    # inlined structure data present (json-escaped into the script)
    assert "receptorPDB" in page and "poseSDF" in page
    from dockstudio import __version__

    assert f"DockStudio v{__version__}" in page


def test_report_neutralizes_script_breakout(tmp_path: Path) -> None:
    receptor_pdb = tmp_path / "rec.pdb"
    receptor_pdb.write_text(
        "REMARK </script><script>alert('XSS')</script>\n"
        "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    pose_sdf = tmp_path / "pose_1.sdf"
    pose_sdf.write_text("test\n  DockStudio\n\n  0  0  0\nM  END\n$$$$\n", encoding="utf-8")
    page = generate_report(
        project=Project(id="p1", name="x", path=str(tmp_path), created_at=datetime.now()),
        request=DockJobRequest(
            receptor_id="r",
            ligand_id="l",
            box=Box(center=(0, 0, 0), size=(1, 1, 1)),
        ),
        result=_sample_result(),
        pose_analyses=[],
        receptor_pdb_path=receptor_pdb,
        pose_sdf_paths=[pose_sdf],
        receptor_name="r",
        ligand_name="l",
    )
    assert "</script><script>alert" not in page
    assert "<\\/script>" in page
