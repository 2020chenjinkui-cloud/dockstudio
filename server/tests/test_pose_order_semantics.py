"""GN-05: poses[0] is engine-order dependent; top_affinity must be order independent.

gnina defaults to --pose_sort_order=CNNscore, so its first pose is the best CNN
score, not the best affinity.  vina writes ascending affinity.  Reading
poses[0].affinity as "the" score therefore compares two different things.
"""
from __future__ import annotations

from dockstudio.models import Box, DockParams, DockResult, Pose


def _result(engine: str, affinities: list[float]) -> DockResult:
    return DockResult(
        engine=engine,
        engine_version="test",
        params=DockParams(),
        box=Box(center=(0, 0, 0), size=(20, 20, 20)),
        poses=[
            Pose(index=i + 1, affinity=a, file=None, rmsd_lb=None, rmsd_ub=None)
            for i, a in enumerate(affinities)
        ],
        run_dir="/tmp/x",
    )


def test_vina_order_is_affinity_order() -> None:
    """vina writes ascending affinity, so poses[0] already is the best."""
    r = _result("vina", [-11.2, -10.4, -9.8])
    assert r.pose_order == "affinity"
    assert r.top_affinity == -11.2
    assert r.poses[0].affinity == r.top_affinity


def test_gnina_order_is_cnn_score_order() -> None:
    """gnina sorts by CNN score: affinity is non-monotonic, poses[0] is NOT the best."""
    r = _result("gnina", [-11.19, -10.47, -11.14])
    assert r.pose_order == "cnn_score"
    assert r.poses[0].affinity == -11.19
    # best affinity is the third pose, not the first
    assert r.top_affinity == -11.19  # -11.19 is still lowest here
    r2 = _result("gnina", [-10.47, -11.19, -11.14])
    assert r2.poses[0].affinity == -10.47
    assert r2.top_affinity == -11.19, "top_affinity must ignore engine order"


def test_top_affinity_ignores_sentinel() -> None:
    """999.0 marks a failed dock and must never be chosen as the top pose."""
    r = _result("gnina", [999.0, -9.5])
    assert r.top_affinity == -9.5


def test_empty_result() -> None:
    r = _result("vina", [])
    assert r.top_pose is None
    assert r.top_affinity is None
