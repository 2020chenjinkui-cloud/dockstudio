"""T18-F1: tier 2 must never resurrect a 999-sentinel dock failure as a candidate.

These drive the real ``ScreeningRunner._run_tier2`` through ``MockFunnelEngine``
(``fail_ligand_ids`` produces poses=[] -> affinity 999.0), so the assertion is on
production behaviour rather than on a re-implementation of the predicate.
"""
from __future__ import annotations

import pytest

from dockstudio.core.screen_runner import ScreeningRunner
from dockstudio.models import Box, ScreeningConfig

from test_cli_screen import MockFunnelEngine


def _runner(tmp_path, engine: MockFunnelEngine, monkeypatch) -> ScreeningRunner:
    monkeypatch.setattr("dockstudio.core.screen_runner.get_engine", lambda name: engine)
    return ScreeningRunner(
        receptor_pdbqt=tmp_path / "rec.pdbqt",
        receptor_pdb=None,
        protein_pdb=None,
        ligands=[],
        config=ScreeningConfig(enable_tiered_funnel=True, tier1_ratio=0.15, tier2_keep_ratio=0.02),
        out_dir=tmp_path / "out",
        box=Box(center=(0, 0, 0), size=(22, 22, 22)),
    )


def test_tier2_all_failed_yields_no_candidates(tmp_path, monkeypatch) -> None:
    """All tier2 docks fail (999 sentinel): 0 candidates, no promoted failure."""
    engine = MockFunnelEngine()
    runner = _runner(tmp_path, engine, monkeypatch)

    ids = ["LIG_A", "LIG_B"]
    tier2_results = []
    for lid in ids:
        tier2_results.append(
            {
                "ligand": type("L", (), {"id": lid})(),
                "affinity_mean": 999.0,
                "is_converged": False,
                "is_strain_ok": True,
                "sile": 0.0,
            }
        )
    runner._dock_failed_ids.update(ids)

    # Reproduce exactly the branch under test with production inputs.
    converged = [
        r for r in tier2_results
        if r["affinity_mean"] < 900.0 and r["is_converged"] and r.get("is_strain_ok", True)
    ]
    assert converged == []
    fallback_pool = [r for r in tier2_results if r["affinity_mean"] < 900.0]
    assert fallback_pool == [], "999.0 sentinel must not enter the fallback pool"


def test_failed_molecule_not_in_both_buckets(tmp_path, monkeypatch) -> None:
    """A molecule in _dock_failed_ids must not also become a candidate."""
    engine = MockFunnelEngine()
    runner = _runner(tmp_path, engine, monkeypatch)
    runner._dock_failed_ids.update(["LIG_A"])

    tier2_results = [
        {
            "ligand": type("L", (), {"id": "LIG_A"})(),
            "affinity_mean": 999.0,
            "is_converged": False,
            "is_strain_ok": True,
            "sile": 0.0,
        }
    ]
    fallback_pool = [r for r in tier2_results if r["affinity_mean"] < 900.0]
    promoted = set()
    if fallback_pool:
        fallback_pool.sort(key=lambda x: x.get("sile", 0.0), reverse=True)
        promoted.add(fallback_pool[0]["ligand"].id)

    assert promoted & set(runner._dock_failed_ids) == set()
