"""gnina/Uni-Dock adapters: availability logic + CLI mapping + output parsing."""

from __future__ import annotations

import platform
import sys
import textwrap
from pathlib import Path

import pytest

from dockstudio.engines.gnina import GninaEngine, parse_gnina_out
from dockstudio.engines.unidock import UniDockEngine
from dockstudio.models import Box, DockParams


def test_gnina_manifest() -> None:
    eng = GninaEngine()
    assert eng.manifest.platforms == ["linux"]
    assert eng.manifest.requires_gpu
    assert eng.manifest.supports_flex_residues
    assert {"vina", "vinardo", "ad4_scoring", "dkoes_fast"} == set(eng.manifest.scoring_functions)


def test_unidock_manifest() -> None:
    eng = UniDockEngine()
    assert eng.manifest.platforms == ["linux"]
    assert eng.manifest.requires_gpu
    assert eng.manifest.supports_batch_ligands
    assert {"vina", "vinardo", "ad4"} == set(eng.manifest.scoring_functions)


@pytest.mark.skipif(platform.system() == "Linux", reason="windows/off-linux check")
def test_gnina_unavailable_reason_chinese() -> None:
    ok, reason = GninaEngine().is_available()
    assert not ok
    assert "WSL2" in (reason or "") and "docs/deploy.md#wsl2" in (reason or "")


@pytest.mark.skipif(platform.system() == "Linux", reason="windows/off-linux check")
def test_unidock_unavailable_reason_chinese() -> None:
    ok, reason = UniDockEngine().is_available()
    assert not ok
    assert "WSL2" in (reason or "")


def test_unidock_cmd_mapping() -> None:
    eng = UniDockEngine()
    box = Box(center=(1, 2, 3), size=(20, 20, 20))
    cmd = eng._build_cmd(
        "unidock", box, DockParams(exhaustiveness=2, num_modes=5), "r.pdbqt", "l.pdbqt", "out.pdbqt"
    )
    assert "--ligands" in cmd and "--ligand" not in cmd
    assert "--exhaustiveness" not in cmd
    i = cmd.index("--search_mode")
    assert cmd[i + 1] == "fast"
    cmd2 = eng._build_cmd("unidock", box, DockParams(exhaustiveness=8), "r", "l", "o")
    assert cmd2[cmd2.index("--search_mode") + 1] == "balance"
    cmd3 = eng._build_cmd("unidock", box, DockParams(exhaustiveness=32), "r", "l", "o")
    assert cmd3[cmd3.index("--search_mode") + 1] == "detail"


def test_gnina_cmd_excludes_energy_range() -> None:
    """gnina CLI 不支持 --energy_range（vina 专有），适配器必须剔除。"""
    from dockstudio.engines.gnina import GninaEngine

    eng = GninaEngine()
    box = Box(center=(1, 2, 3), size=(20, 20, 20))
    cmd = eng._build_cmd(
        "gnina", box, DockParams(exhaustiveness=8, energy_range=3.0), "r", "l", "o"
    )
    assert "--energy_range" not in cmd
    assert cmd[cmd.index("--exhaustiveness") + 1] == "8"


def test_gnina_cpu_opt_in_flag(monkeypatch) -> None:
    from dockstudio.engines.gnina import GninaEngine

    monkeypatch.setenv("DOCKSTUDIO_GNINA_CPU", "1")
    eng = GninaEngine()
    cmd = eng._build_cmd(
        "gnina",
        Box(center=(0, 0, 0), size=(20, 20, 20)),
        DockParams(exhaustiveness=1),
        "r",
        "l",
        "o",
    )
    assert "--no_gpu" in cmd


GNINA_OUT = textwrap.dedent(
    """\
    MODEL 1
    REMARK minimizedAffinity -9.5
    ROOT
    ATOM      1  C   UNL     1       1.0 2.0 3.0
    ENDROOT
    ENDMDL
    MODEL 2
    REMARK minimizedAffinity -8.7
    ROOT
    ATOM      1  C   UNL     1       2.0 2.0 3.0
    ENDROOT
    ENDMDL
    """
)

GNINA_LOG = (
    "Using random seed: 42\nCNNscore 0.98 CNNaffinity -9.9\nCNNscore 0.91 CNNaffinity -8.9\n"
)


def test_gnina_parse_remarks(tmp_path: Path) -> None:
    out = tmp_path / "out.pdbqt"
    out.write_text(GNINA_OUT, encoding="utf-8")
    poses = parse_gnina_out(out, GNINA_LOG)
    assert len(poses) == 2
    assert poses[0].affinity == pytest.approx(-9.5)
    assert poses[1].affinity == pytest.approx(-8.7)


def test_gnina_parse_garbage_raises(tmp_path: Path) -> None:
    from dockstudio.models import DockingError

    out = tmp_path / "out.pdbqt"
    out.write_text("junk\n", encoding="utf-8")
    with pytest.raises(DockingError):
        parse_gnina_out(out, "nothing")


@pytest.mark.skipif(
    sys.platform == "win32",
    reason=(
        "T12-F1: fake-gnina stub is a #!<python> script; Windows CreateProcess "
        "cannot execute it (no shell, shebang ignored). Covered on Linux/macOS."
    ),
)
def test_gnina_mock_binary_e2e(tmp_path: Path, monkeypatch) -> None:
    """Whole dock() flow with a fake gnina shell script binary."""
    # T9-F1: a Python stub with an "env python" shebang instead of a /bin/sh
    # script - Windows CI has no `sh`, and the old stub made this file fail on
    # the windows-latest leg of server-test-full.
    script = tmp_path / "gnina"
    script.write_text(
        f"#!{sys.executable}\n"
        "import sys\n"
        "from pathlib import Path\n"
        f"OUT = {GNINA_OUT!r}\n"
        "args = sys.argv[1:]\n"
        "if '--version' in args:\n"
        "    print('gnina 1.0 fake')\n"
        "    raise SystemExit(0)\n"
        "prev = ''\n"
        "out = None\n"
        "for a in args:\n"
        "    if prev == '--out':\n"
        "        out = a\n"
        "    prev = a\n"
        "if out:\n"
        "    Path(out).write_text(OUT, encoding='utf-8')\n"
        "print('1  -9.5  0.0  0.0')\n"
        "print('CNNscore 0.98 CNNaffinity -9.9')\n"
        "raise SystemExit(0)\n",
        encoding="utf-8",
    )
    script.chmod(0o755)
    monkeypatch.setenv("DOCKSTUDIO_GNINA", str(script))
    monkeypatch.setattr(platform, "system", lambda: "Linux")
    monkeypatch.setattr(GninaEngine, "_gpu_ok", lambda self: (True, ""))
    monkeypatch.setattr(GninaEngine, "_version", lambda self, binary: "1.0-mock")

    eng = GninaEngine()
    ok, reason = eng.is_available()
    assert ok, reason
    result = None
    try:
        work = tmp_path / "run"
        rec = tmp_path / "rec.pdbqt"
        lig = tmp_path / "lig.pdbqt"
        # The receptor must look like a real PDBQT: atom-less input is now
        # rejected by the adapter itself (GN-04), and this test is about the
        # command/parse flow, not about input validation.
        rec.write_text(
            "ATOM      1  N   ILE A  10      13.000  22.000   5.000  1.00  0.00    -0.300 N\n"
            "ATOM      2  CA  ILE A  10      14.000  22.500   5.200  1.00  0.00     0.100 C\n",
            encoding="utf-8",
        )
        lig.write_text("x", encoding="utf-8")
        result = eng.dock(
            rec,
            lig,
            Box(center=(0, 0, 0), size=(20, 20, 20)),
            DockParams(exhaustiveness=2, num_modes=2, seed=42),
            work,
        )
    except OSError as exc:
        # T15-F4: only a *launch* failure (ENOEXEC / EACCES / ENOENT) means "this
        # host cannot run a script stub" - that is a skip.  Anything else is a real
        # failure of the code under test and must surface as such, not be relabelled
        # as "not applicable here" (which is the same disease as the swallowed
        # aborts this round is about).
        pytest.skip(f"host cannot launch a script stub: {exc}")
    assert result is not None
    assert result.poses[0].affinity == pytest.approx(-9.5)
    assert "CNNscore" in result.log


def test_unidock_build_batch_cmd() -> None:
    """Verifies UniDockEngine._build_batch_cmd builds exact --ligand_index and --dir CLI."""
    eng = UniDockEngine()
    box = Box(center=(10.5, 20.0, -5.2), size=(22.0, 24.0, 20.0))
    params = DockParams(exhaustiveness=4, num_modes=9, scoring="vina", seed=2026)

    cmd = eng._build_batch_cmd(
        binary="unidock",
        box=box,
        params=params,
        receptor_path="/path/rec.pdbqt",
        ligand_index_path="/path/index.txt",
        out_dir_path="/path/out_dir",
    )

    assert "--receptor" in cmd and cmd[cmd.index("--receptor") + 1] == "/path/rec.pdbqt"
    assert "--ligand_index" in cmd and cmd[cmd.index("--ligand_index") + 1] == "/path/index.txt"
    assert "--dir" in cmd and cmd[cmd.index("--dir") + 1] == "/path/out_dir"
    assert "--search_mode" in cmd and cmd[cmd.index("--search_mode") + 1] == "fast"
    assert "--num_modes" in cmd and cmd[cmd.index("--num_modes") + 1] == "9"
    assert "--scoring" in cmd and cmd[cmd.index("--scoring") + 1] == "vina"
    assert "--seed" in cmd and cmd[cmd.index("--seed") + 1] == "2026"

    # Assert single-ligand / inapplicable flags are excluded
    assert "--ligand" not in cmd
    assert "--ligands" not in cmd
    assert "--out" not in cmd
    assert "--cpu" not in cmd
    assert "--energy_range" not in cmd
