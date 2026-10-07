"""Geometric pocket detection on a 1 Å grid (pure numpy/scipy, no OpenBabel).

Method ("builtin_grid"): rasterize protein heavy atoms into a grid padded by
4 Å; free voxels are those 2.0–6.5 Å from the nearest atom (clear of the
protein but still near its surface). Each free voxel gets a *burial index* =
number of protein atoms within a 12 Å sphere (a depth/concavity proxy — high
inside pocket mouths and buried grooves, low on flat surfaces). Voxels above
the 90th percentile are clustered by 6-connectivity; clusters ≥ 20 voxels are
reported as pockets ranked by voxel count. Each pocket box = cluster envelope
with 4 Å padding, clamped to 15–40 Å per side.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from ..logging import get_logger
from ..models import Pocket, PocketResponse, PreparationError

log = get_logger("analysis.pockets")

_GRID_STEP = 1.0
_CLEAR_MIN = 2.0  # voxels closer than this to an atom are protein volume
_CLEAR_MAX = 6.5  # voxels farther than this are solvent, not pocket candidates
_BURIAL_R = 12.0  # sphere radius for the burial index
_BURIAL_Q = 0.90  # keep free voxels above this burial-index quantile
_MIN_POCKET_POINTS = 20
_BOX_MIN = 15.0
_BOX_MAX = 40.0
_PADDING = 4.0
_MARGIN = 4.0


def _load_protein_coords(pdb_path: Path) -> np.ndarray:
    coords: list[list[float]] = []
    for line in pdb_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("ATOM  ") and len(line) >= 54:
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            if not element:
                element = line[12:16].strip().lstrip("0123456789")[:1].upper()
            if element == "H":
                continue
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
    if not coords:
        raise PreparationError(f"no protein atoms found in {pdb_path.name}")
    return np.asarray(coords, dtype=np.float64)


def detect_pockets(pdb_path: Path, top_n: int = 5) -> list[Pocket]:
    """Detect pockets in a protein PDB; returns ranked pocket boxes."""
    coords = _load_protein_coords(pdb_path)
    tree = cKDTree(coords)

    lo = coords.min(axis=0) - _MARGIN
    hi = coords.max(axis=0) + _MARGIN
    shape = np.ceil((hi - lo) / _GRID_STEP).astype(int) + 1
    atom_grid = np.zeros(shape, dtype=bool)
    idx = np.clip(np.round((coords - lo) / _GRID_STEP).astype(int), 0, shape - 1)
    atom_grid[tuple(idx.T)] = True
    # distance (in Å) from each grid voxel to the nearest protein atom
    dist = ndimage.distance_transform_edt(~atom_grid) * _GRID_STEP

    free_idx = np.argwhere((dist >= _CLEAR_MIN) & (dist < _CLEAR_MAX))
    if len(free_idx) == 0:
        return []
    world = lo + free_idx * _GRID_STEP
    burial = tree.query_ball_point(world, r=_BURIAL_R, return_length=True)
    threshold = max(float(np.quantile(burial, _BURIAL_Q)), 80.0)
    candidates = free_idx[burial >= threshold]

    mask = np.zeros(shape, dtype=bool)
    if len(candidates):
        mask[tuple(candidates.T)] = True
    structure6 = ndimage.generate_binary_structure(3, 1)
    labels, n_clusters = ndimage.label(mask, structure=structure6)

    pockets: list[Pocket] = []
    for cid in range(1, n_clusters + 1):
        member = np.argwhere(labels == cid)
        if len(member) < _MIN_POCKET_POINTS:
            continue
        w = lo + member * _GRID_STEP
        center = w.mean(axis=0)
        size = np.clip((w.max(axis=0) + _PADDING) - (w.min(axis=0) - _PADDING), _BOX_MIN, _BOX_MAX)
        pockets.append(
            Pocket(
                id=0,
                center=tuple(round(float(c), 3) for c in center),  # type: ignore[arg-type]
                size=tuple(round(float(s), 3) for s in size),  # type: ignore[arg-type]
                score=0.0,
                n_points=len(member),
            )
        )
    pockets.sort(key=lambda p: p.n_points, reverse=True)
    top = pockets[: max(1, top_n)]
    if top:
        max_pts = top[0].n_points
        for rank, p in enumerate(top, start=1):
            p.id = rank
            p.score = round(p.n_points / max_pts, 3)
    log.info("pockets in %s: %d kept of %d clusters", pdb_path.name, len(top), n_clusters)
    return top


def pocket_response(
    receptor_id: str, protein_pdb: Path, cache_dir: Path, top_n: int = 5
) -> PocketResponse:
    """Detect pockets with per-(receptor_id, top_n) JSON caching."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"pockets_top{top_n}.json"
    cached = _load_cache(cache_path, receptor_id)
    if cached is not None:
        return cached
    pockets = detect_pockets(protein_pdb, top_n=top_n)
    response = PocketResponse(receptor_id=receptor_id, pockets=pockets)
    cache_path.write_text(response.model_dump_json(indent=2), encoding="utf-8")
    return response


def _load_cache(cache_path: Path, receptor_id: str) -> PocketResponse | None:
    if not cache_path.is_file():
        return None
    try:
        data = json.loads(cache_path.read_text(encoding="utf-8"))
        if data.get("receptor_id") == receptor_id:
            return PocketResponse(**data)
    except (OSError, json.JSONDecodeError, ValueError):
        pass
    return None
