"""打包回归：PyInstaller spec 必须收集 demo_data，否则冻结包里 /api/demo/example 直接 FileNotFoundError。

生产事故复盘：onefile 桌面后端缺 1HSG 演示数据 → 500 → 前端报"后端未启动"。
"""

from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_onefile_spec_collects_demo_data() -> None:
    spec = (REPO_ROOT / "scripts" / "backend_onefile.spec").read_text(encoding="utf-8")
    assert "demo_data" in spec, "backend_onefile.spec 未收集 dockstudio/demo_data"
    assert "plip_worker.py" in spec, "backend_onefile.spec 未收集 PLIP worker"
    assert "meeko.cli.mk_prepare_receptor" in spec, "缺少冻结端可用的 Meeko CLI"
    assert "pyi_rdkit_data_hook.py" in spec, "缺少非 ASCII 路径的 RDKit 数据 hook"
    assert 'collect_data_files("posebusters")' in spec, "冻结端缺少 PoseBusters 配置 yml"


def test_folder_spec_collects_demo_data() -> None:
    spec = (REPO_ROOT / "scripts" / "backend.spec").read_text(encoding="utf-8")
    assert "demo_data" in spec, "backend.spec 未收集 dockstudio/demo_data"
    assert "plip_worker.py" in spec, "backend.spec 未收集 PLIP worker"
    assert "meeko.cli.mk_prepare_receptor" in spec, "缺少冻结端可用的 Meeko CLI"
    assert "pyi_rdkit_data_hook.py" in spec, "缺少非 ASCII 路径的 RDKit 数据 hook"
    assert 'collect_data_files("posebusters")' in spec, "冻结端缺少 PoseBusters 配置 yml"


def test_demo_data_files_present() -> None:
    demo = REPO_ROOT / "server" / "dockstudio" / "demo_data"
    assert (demo / "1HSG.pdb").is_file()
    assert (demo / "indinavir_smiles.txt").is_file()
