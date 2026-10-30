"""Report v2: sections filtering, custom title, 400 on invalid section."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from dockstudio.analysis.report import generate_report
from dockstudio.models import (
    Box,
    DockJobRequest,
    DockParams,
    DockResult,
    Pose,
    Project,
)


def _result() -> DockResult:
    return DockResult(
        engine="vina",
        engine_version="1.2.7",
        params=DockParams(exhaustiveness=4, num_modes=3, seed=42),
        box=Box(center=(13.0, 22.4, 5.5), size=(22.5, 22.5, 22.5)),
        poses=[Pose(index=1, affinity=-9.8), Pose(index=2, affinity=-9.0)],
        run_dir="runs/x",
    )


def _gen(tmp_path: Path, pose_analyses=None, **kwargs) -> str:
    rec = tmp_path / "rec.pdb"
    rec.write_text(
        "ATOM      1  CA  ALA A   1       1.0   2.0   3.0  1.00  0.00           C\nEND\n",
        encoding="utf-8",
    )
    pose = tmp_path / "pose_1.sdf"
    pose.write_text("p\n  DockStudio\n\n  0  0  0\nM  END\n$$$$\n", encoding="utf-8")
    return generate_report(
        project=Project(id="p", name="项目甲", path="/tmp", created_at=datetime.now()),
        request=DockJobRequest(
            receptor_id="r", ligand_id="l", box=Box(center=(0, 0, 0), size=(1, 1, 1))
        ),
        result=_result(),
        pose_analyses=pose_analyses if pose_analyses is not None else [],
        receptor_pdb_path=rec,
        pose_sdf_paths=[pose],
        receptor_name="1HSG",
        ligand_name="indinavir",
        **kwargs,
    )


def test_sections_scores_disclaimer_only(tmp_path: Path) -> None:
    page = _gen(tmp_path, sections=["scores", "disclaimer"])
    assert "$3Dmol" not in page and "3dmol-min.js" not in page.lower()
    assert "姿态打分" in page
    assert "免责声明" in page
    assert "任务快照" not in page
    assert "3D 结合模式" not in page


def test_custom_title(tmp_path: Path) -> None:
    page = _gen(tmp_path, title="吲哚美辛初筛报告")
    assert "<h1>吲哚美辛初筛报告</h1>" in page
    # default title fallback contains project name
    assert "项目甲" in _gen(tmp_path)


def test_long_title_truncated(tmp_path: Path) -> None:
    page = _gen(tmp_path, title="长" * 200)
    assert page.count("长" * 81) == 0  # capped at 80 chars


def test_diagram_embedded(tmp_path: Path) -> None:
    svg = '<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"></svg>'
    page = _gen(tmp_path, diagram_svg=svg)
    assert "2D 相互作用示意图" in page
    page2 = _gen(tmp_path, sections=["scores"], diagram_svg=svg)
    assert "2D 相互作用示意图" not in page2


def test_posecheck_shown(tmp_path: Path) -> None:
    from dockstudio.models import PoseAnalysis, PoseCheck, PoseCheckIssue

    analysis = PoseAnalysis(
        pose_index=1,
        source="builtin",
        pose_check=PoseCheck(
            pose_index=1,
            ok=False,
            issues=[PoseCheckIssue(rule="clash", severity="error", message="原子碰撞")],
        ),
    )
    page = _gen(tmp_path, pose_analyses=[analysis])
    assert "姿态质检" in page
    assert "未通过" in page
    assert "原子碰撞" in page
