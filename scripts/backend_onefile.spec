# -*- mode: python ; coding: utf-8 -*-
"""onefile 变体：单一可执行 dockstudio-backend(.exe)，供 Tauri sidecar 携带。
与 backend.spec 同源的收集规则；onefile 启动时解压较慢（科学计算库大，5-15s），
换来的是桌面安装包只带一个后端文件。
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files

ROOT = Path(SPECPATH).parent  # noqa: F821

block_cipher = None

_ob_datas, _ob_binaries, _ob_hiddenimports = collect_all("openbabel")
_rdkit_datas = collect_data_files("rdkit")
_meeko_datas = collect_data_files("meeko")
_gemmi_datas = collect_data_files("gemmi")

import importlib.util as _ilu

_ob_spec = _ilu.find_spec("openbabel")
_ob_pkg = Path(_ob_spec.submodule_search_locations[0]) if _ob_spec else None
_ob_explicit_datas: list[tuple[str, str]] = []
if _ob_pkg is not None:
    for _sub in ("bin", "lib", "share"):
        if (_ob_pkg / _sub).is_dir():
            _ob_explicit_datas.append((str(_ob_pkg / _sub), f"openbabel/{_sub}"))

a = Analysis(
    [str(ROOT / "scripts" / "server_entry.py")],
    pathex=[str(ROOT / "server")],
    binaries=_ob_binaries,
    datas=[
        (str(ROOT / "resources" / "bin"), "resources/bin"),
        # 内置演示项目数据（1HSG 受体/配体/SMILES），冻结包里缺了 /api/demo/example 直接 500
        (str(ROOT / "server" / "dockstudio" / "demo_data"), "dockstudio/demo_data"),
    ]
    + _ob_datas
    + _rdkit_datas
    + _meeko_datas
    + _gemmi_datas
    + _ob_explicit_datas,
    hiddenimports=[
        "rdkit",
        "rdkit.Chem.AllChem",
        "meeko",
        "openbabel",
        "scipy.spatial",
        "gemmi",
        "uvicorn.lifespan",
        "uvicorn.protocols",
    ]
    + _ob_hiddenimports,
    hookspath=[],
    runtime_hooks=[str(ROOT / "scripts" / "pyi_obabel_hook.py")],
    excludes=["tkinter", "matplotlib"],
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name="dockstudio-backend",
    debug=False,
    strip=False,
    upx=False,
    console=False,
    icon=str(ROOT / "apps" / "desktop" / "src-tauri" / "icons" / "icon.ico"),
)
