# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec：把 FastAPI 后端 + 科学计算依赖 + Vina 二进制打成单目录产物。

用法（各平台 CI 上执行）：
    py -3.12 -m pip install pyinstaller          # Windows；Linux/mac 用 python3 -m pip
    py -3.12 -m PyInstaller scripts/backend.spec --noconfirm

产物 dist/dockstudio-backend/ 内含可执行 dockstudio-backend(.exe)，
Tauri 将其作为 sidecar 携带；也可独立分发作免前端的服务器模式。
注意：openbabel 的格式插件（BABEL_LIBDIR）与 rdkit 的数据文件（BaseFeatures 等）
必须随包收集，否则冻结后端启动即崩（macOS CI 实测）。
"""

from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files

ROOT = Path(SPECPATH).parent  # noqa: F821  SPECPATH 由 PyInstaller 注入

block_cipher = None

_ob_datas, _ob_binaries, _ob_hiddenimports = collect_all("openbabel")
_rdkit_datas = collect_data_files("rdkit")
_meeko_datas = collect_data_files("meeko")  # ad4_types 等参数 JSON，缺了配体准备即崩
_gemmi_datas = collect_data_files("gemmi")

# openbabel 插件/数据目录强制随包（collect_all 对平台特定子布局不保险，
# mac CI 实测 bin/lib/share 收集不全导致 LoadAllPlugins 失败）
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
    [],
    exclude_binaries=True,
    name="dockstudio-backend",
    debug=False,
    strip=False,
    upx=False,
    console=True,
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    name="dockstudio-backend",
)
