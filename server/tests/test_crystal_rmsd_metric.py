"""D9 metric semantics: in-place (CalcRMS) vs superposed (GetBestRMS).

D9 (交付物清单「Redocking 对照」) requires the co-crystal ligand to be re-docked
*at its original site*, i.e. a pose that has drifted out of the pocket must be
reported as a failure.  A metric that first superimposes the pose onto the
reference cannot express that: it measures conformation, not placement.

This module pins the contract with a geometry where the two answers differ by
a whole order of magnitude: the same molecule twice, the `pose` translated by
+10 A in x relative to the `crystal`.  Internal geometry is identical, so

    superposed / symmetry-corrected  -> ~0.0 A   (GetBestRMS)
    in-place   / symmetry-corrected  -> ~10.0 A  (CalcRMS)

The pair is embedded inline (produced once with ETKDGv3/MMFF94s, randomSeed=42,
phenylacetic acid C8H8O2, 10 heavy atoms) so the test has no embedder dependency.
"""

from __future__ import annotations

from pathlib import Path

import pytest

# --- synthetic pair: same molecule, pose translated +10 A along x -------------
_CRYSTAL_BLOCK = """phenylacetic acid
     RDKit          3D

 18 18  0  0  0  0  0  0  0  0999 V2000
   -0.5679    1.2793   -0.1002 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8929    1.2984    0.3377 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.6108    0.1086    0.4429 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.0036   -1.1020    0.1143 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6784   -1.1237   -0.3233 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0493    0.0682   -0.4422 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4847    0.0450   -0.8970 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.4563    0.0497    0.2520 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8492    1.0291    0.8660 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.8784   -1.1850    0.5935 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0173    2.2158   -0.1670 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.3636    2.2422    0.6006 H   0  0  0  0  0  0  0  0  0  0  0  0
   -3.6421    0.1245    0.7854 H   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5615   -2.0306    0.2019 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2165   -2.0779   -0.5683 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6999    0.9138   -1.5310 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6780   -0.8279   -1.5327 H   0  0  0  0  0  0  0  0  0  0  0  0
    3.4589   -1.0276    1.3675 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  7  1  0
  7  8  1  0
  8  9  2  0
  8 10  1  0
  6  1  2  0
  1 11  1  0
  2 12  1  0
  3 13  1  0
  4 14  1  0
  5 15  1  0
  7 16  1  0
  7 17  1  0
 10 18  1  0
M  END
"""


def _counts_line_index(lines: list[str]) -> int:
    """Locate the V2000 counts line (title/program/comment precede it)."""
    for i, line in enumerate(lines[:6]):
        if "V2000" in line and line[0:3].strip().isdigit():
            return i
    raise AssertionError("not a V2000 molblock")


def _translate(block: str, dx: float) -> str:
    """Shift every atom's x coordinate by dx, keeping the molblock otherwise intact.

    V2000: the counts line is followed by one fixed-width line per atom
    (`%10.4f%10.4f%10.4f` for x/y/z).  The header length varies because RDKit
    writes an empty title line, so the counts line is located rather than assumed.
    """
    lines = block.splitlines()
    counts = _counts_line_index(lines)
    natoms = int(lines[counts][0:3])
    shifted = list(lines)
    for i in range(counts + 1, counts + 1 + natoms):
        x = float(lines[i][0:10]) + dx
        shifted[i] = f"{x:10.4f}" + lines[i][10:]
    return "\n".join(shifted)


@pytest.fixture()
def pair(tmp_path: Path) -> tuple[Path, Path]:
    block = _CRYSTAL_BLOCK.strip("\n")
    crystal = tmp_path / "crystal.sdf"
    crystal.write_text(block + "\n$$$$\n", encoding="utf-8")
    pose = tmp_path / "pose_translated.sdf"
    pose.write_text(_translate(block, 10.0) + "\n$$$$\n", encoding="utf-8")
    return crystal, pose


# ---------------------------------------------------------------------------
# the contract: two metrics, two different questions, two different answers
# ---------------------------------------------------------------------------
def test_translated_pose_inplace_is_ten_angstrom(pair) -> None:
    """D9 judges placement: a pose 10 A away must score ~10 A, not ~0 A."""
    from dockstudio.analysis.crystal import crystal_pose_rmsd

    crystal, pose = pair
    result = crystal_pose_rmsd(crystal, pose)
    assert result.inplace_rmsd is not None
    assert 9.5 <= result.inplace_rmsd <= 10.5, (
        f"in-place RMSD should be ~10 A for a pose translated by 10 A, "
        f"got {result.inplace_rmsd}"
    )


def test_translated_pose_superposed_is_near_zero(pair) -> None:
    """The superposed metric is the *auxiliary* one: it only sees conformation."""
    from dockstudio.analysis.crystal import crystal_pose_rmsd

    crystal, pose = pair
    result = crystal_pose_rmsd(crystal, pose)
    assert result.superposed_rmsd is not None
    assert result.superposed_rmsd < 0.5, (
        f"superposed RMSD should stay near 0 for a pure translation, "
        f"got {result.superposed_rmsd}"
    )


def test_identical_pose_scores_zero_on_both(pair) -> None:
    """Sanity: no drift means both metrics agree near zero."""
    from dockstudio.analysis.crystal import crystal_pose_rmsd

    crystal, _pose = pair
    result = crystal_pose_rmsd(crystal, crystal)
    assert result.inplace_rmsd is not None and result.superposed_rmsd is not None
    assert result.inplace_rmsd < 0.05
    assert result.superposed_rmsd < 0.05


def test_metric_names_are_distinguishable(pair) -> None:
    """Callers must be able to tell which number is the D9 verdict."""
    from dockstudio.analysis.crystal import (
        CRITERION_METRIC,
        INPLACE_METRIC,
        SUPERPOSED_METRIC,
        crystal_pose_rmsd,
    )

    assert CRITERION_METRIC == INPLACE_METRIC
    assert CRITERION_METRIC != SUPERPOSED_METRIC
    crystal, pose = pair
    result = crystal_pose_rmsd(crystal, pose)
    assert result.metric == CRITERION_METRIC
    assert result.auxiliary_metric == SUPERPOSED_METRIC
    # the D9 number is the in-place one, never the superposed one
    assert result.rmsd == result.inplace_rmsd
    assert result.rmsd != result.superposed_rmsd


def test_d9_pass_uses_inplace_metric(pair) -> None:
    """A 10 A translation must fail D9 even though the geometry is identical."""
    from dockstudio.analysis.crystal import D9_THRESHOLD_A, crystal_pose_rmsd

    crystal, pose = pair
    far = crystal_pose_rmsd(crystal, pose)
    near = crystal_pose_rmsd(crystal, crystal)
    assert far.d9_pass is False, "in-place RMSD ~10 A must not pass D9"
    assert near.d9_pass is True
    assert far.threshold_a == D9_THRESHOLD_A == 2.0


def test_heavy_atom_rmsd_docstring_states_its_metric() -> None:
    """The legacy entry point keeps working but must declare which metric it is."""
    from dockstudio.analysis import crystal as crystal_mod

    doc = (crystal_mod.heavy_atom_rmsd.__doc__ or "").lower()
    assert "superimpos" in doc or "superposed" in doc or "叠合" in doc