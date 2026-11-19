"""D9 wiring: the in-place verdict must reach the API, the report and the audit bundle.

The metric-level contract lives in ``test_crystal_rmsd_metric.py``; this module
checks the *consumers* — a reader must always be able to tell which function
produced a number, and the D9 verdict must be the in-place one.
"""

from __future__ import annotations

import json
import time
import zipfile
from datetime import datetime
from pathlib import Path

import pytest

from dockstudio.analysis.crystal import (
    CRITERION_METRIC,
    D9_THRESHOLD_A,
    INPLACE_METRIC,
    SUPERPOSED_METRIC,
    crystal_pose_rmsd,
)
from dockstudio.analysis.report import ALL_SECTIONS, generate_report
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


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
def _gen(tmp_path: Path, **kwargs) -> str:
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
        pose_analyses=[],
        receptor_pdb_path=rec,
        pose_sdf_paths=[pose],
        receptor_name="1HSG",
        ligand_name="indinavir",
        **kwargs,
    )


def test_redocking_section_is_registered() -> None:
    assert "redocking" in ALL_SECTIONS


def test_report_shows_both_metrics_with_names(tmp_path: Path) -> None:
    from dockstudio.analysis.crystal import CrystalPoseRmsd

    series = {
        1: CrystalPoseRmsd(inplace_rmsd=0.454, superposed_rmsd=0.380),
        2: CrystalPoseRmsd(inplace_rmsd=10.9, superposed_rmsd=0.43),
    }
    page = _gen(tmp_path, crystal_reference="MK1:B", crystal_rmsd_series=series)
    assert "共晶配体自对接对照（D9）" in page
    assert "0.454" in page and "0.380" in page
    assert "10.900" in page
    # verdicts: pose 1 passes, pose 2 fails — decided by the in-place number
    assert "PASS" in page and "FAIL" in page
    # both metric identifiers are spelled out so the numbers cannot be confused
    assert INPLACE_METRIC in page
    assert SUPERPOSED_METRIC in page
    assert "不叠合" in page


def test_report_best_of_n_uses_inplace(tmp_path: Path) -> None:
    """Best-of-N must be picked by the D9 metric, not by the superposed one."""
    from dockstudio.analysis.crystal import CrystalPoseRmsd

    series = {
        1: CrystalPoseRmsd(inplace_rmsd=4.5, superposed_rmsd=0.10),  # great shape, wrong place
        2: CrystalPoseRmsd(inplace_rmsd=0.60, superposed_rmsd=0.55),
    }
    page = _gen(tmp_path, crystal_reference="MK1:B", crystal_rmsd_series=series)
    assert "best-of-2" in page
    assert "0.600" in page, "best-of-N must report the best in-place value"


def test_report_redocking_absent_without_series(tmp_path: Path) -> None:
    page = _gen(tmp_path)
    assert "该受体没有可用的共晶配体" in page
    assert "SKIP" in page


def test_report_redocking_section_can_be_filtered(tmp_path: Path) -> None:
    from dockstudio.analysis.crystal import CrystalPoseRmsd

    series = {1: CrystalPoseRmsd(inplace_rmsd=0.5, superposed_rmsd=0.4)}
    page = _gen(tmp_path, sections=["scores"], crystal_reference="MK1:B", crystal_rmsd_series=series)
    assert "共晶配体自对接对照" not in page


# ---------------------------------------------------------------------------
# API endpoint
# ---------------------------------------------------------------------------
def test_api_crystal_rmsd_reports_both_metrics(client) -> None:
    """The endpoint exposes the in-place verdict plus the auxiliary metric."""
    from test_audit import _seed  # reuse the fake-engine fixture helper

    monkeypatch = pytest.MonkeyPatch()
    try:
        jid, _pid, rec, _lig = _seed(client, monkeypatch)
    finally:
        monkeypatch.undo()
    resp = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 1})
    assert resp.status_code == 200
    data = resp.json()
    assert data["reference"] == "MK1:B"
    # both numbers present and distinguishable
    assert "rmsd" in data and "rmsd_superposed" in data
    assert data["metric"] == CRITERION_METRIC == INPLACE_METRIC
    assert data["metric_superposed"] == SUPERPOSED_METRIC
    assert data["threshold_a"] == D9_THRESHOLD_A
    assert data["d9_status"] in {"PASS", "FAIL", "SKIP"}
    # the D9 verdict follows the in-place value, never the superposed one
    if data["rmsd"] is not None:
        assert data["d9_status"] == ("PASS" if data["rmsd"] <= D9_THRESHOLD_A else "FAIL")


def test_api_crystal_rmsd_null_case_keeps_metric_names(client) -> None:
    """No co-crystal: numbers are null but the metric names still identify the口径."""
    from test_audit import _seed

    monkeypatch = pytest.MonkeyPatch()
    try:
        jid, _pid, rec, _lig = _seed(client, monkeypatch)
    finally:
        monkeypatch.undo()
    meta = Path(rec.pdbqt_path).parent / "cocrystals.json"
    meta.write_text("[]", encoding="utf-8")
    resp = client.get(f"/api/jobs/{jid}/crystal_rmsd", params={"pose_index": 1})
    data = resp.json()
    assert data["rmsd"] is None and data["rmsd_superposed"] is None
    assert data["metric"] == INPLACE_METRIC
    assert data["metric_superposed"] == SUPERPOSED_METRIC
    assert data["d9_status"] == "SKIP"


# ---------------------------------------------------------------------------
# audit bundle
# ---------------------------------------------------------------------------
def test_audit_bundle_carries_d9_section(client) -> None:
    from test_audit import _seed

    monkeypatch = pytest.MonkeyPatch()
    try:
        jid, _pid, _rec, _lig = _seed(client, monkeypatch)
    finally:
        monkeypatch.undo()
    resp = client.get(f"/api/jobs/{jid}/audit.zip")
    assert resp.status_code == 200
    tmp = Path(client.app.state.settings.data_dir) / "d9_audit.zip"
    tmp.write_bytes(resp.content)
    with zipfile.ZipFile(tmp) as zf:
        report = zf.read("report.md").decode("utf-8")
        snapshot = json.loads(zf.read("snapshot.json"))
    assert "共晶配体自对接对照 (D9)" in report
    assert INPLACE_METRIC in report and SUPERPOSED_METRIC in report
    assert "不叠合" in report
    assert "不可作为 D9 结论" in report
    # snapshot.json keeps the machine-readable form
    assert "cocrystal_rmsd" in snapshot
    block = snapshot["cocrystal_rmsd"]
    assert block["criterion_metric"] == INPLACE_METRIC
    assert block["auxiliary_metric"] == SUPERPOSED_METRIC
    assert block["threshold_a"] == D9_THRESHOLD_A
    for pose_payload in block["poses"].values():
        assert pose_payload["metric"] == INPLACE_METRIC
        assert pose_payload["rmsd"] == pose_payload["rmsd"]
        assert pose_payload["auxiliary_metric"] == SUPERPOSED_METRIC


# ---------------------------------------------------------------------------
# end-to-end: a real docked pose must be judged by placement
# ---------------------------------------------------------------------------
def _write_sdf_pair(tmp_path: Path) -> tuple[Path, Path]:
    from test_crystal_rmsd_metric import _CRYSTAL_BLOCK, _translate

    block = _CRYSTAL_BLOCK.strip("\n")
    crystal = tmp_path / "crystal.sdf"
    crystal.write_text(block + "\n$$$$\n", encoding="utf-8")
    far = tmp_path / "pose_far.sdf"
    far.write_text(_translate(block, 10.0) + "\n$$$$\n", encoding="utf-8")
    return crystal, far


def test_series_and_best_of_n_prefer_placement(tmp_path: Path) -> None:
    """A well-shaped but off-site pose must not win best-of-N."""
    from dockstudio.analysis.crystal import best_pose_rmsd, pose_rmsd_series

    crystal, far = _write_sdf_pair(tmp_path)
    near = tmp_path / "pose_near.sdf"
    near.write_text(crystal.read_text(encoding="utf-8"), encoding="utf-8")

    series = pose_rmsd_series(crystal, [far, near])
    assert series[1].inplace_rmsd > 9.0  # off-site
    assert series[1].superposed_rmsd < 0.5  # ... but the same conformation
    assert best_pose_rmsd(series) is series[2]


def test_cocrystal_reference_reads_sidecar(tmp_path: Path) -> None:
    from dockstudio.analysis.crystal import cocrystal_reference

    assert cocrystal_reference(tmp_path) == (None, None)
    sdf = tmp_path / "cocrystal_MK1_B.sdf"
    sdf.write_text("x\n  DockStudio\n\n  0  0  0\nM  END\n$$$$\n", encoding="utf-8")
    (tmp_path / "cocrystals.json").write_text(
        json.dumps([{"resname": "MK1", "chain": "B", "sdf_path": str(sdf)}]), encoding="utf-8"
    )
    reference, path = cocrystal_reference(tmp_path)
    assert reference == "MK1:B"
    assert path == sdf
    # a missing file degrades to "no reference" rather than raising
    sdf.unlink()
    assert cocrystal_reference(tmp_path) == ("MK1:B", None)