# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec：把 FastAPI 后端 + 科学计算依赖 + Vina 二进制打成单目录产物。

用法（各平台 CI 上执行）：
    py -3.12 -m pip install pyinstaller          # Windows；Linux/mac 用 python3 -m pip
    py -3.12 -m PyInstaller scripts/backend.spec --noconfirm

产物 dist/dockstudio-backend/ 内含可执行 dockstudio-backend(.exe)，
Tauri 将其作为 sidecar 携带；也可独立分发作免前端的服务器模式。
注意：rdkit/openbabel/meeko 的隐式 import 若漏打包，在 hiddenimports 里补。
"""

from pathlib import Path

ROOT = Path(SPECPATH).parent  # noqa: F821  SPECPATH 由 PyInstaller 注入

block_cipher = None

a = Analysis(
    [str(ROOT / "scripts" / "server_entry.py")],
    pathex=[str(ROOT / "server")],
    binaries=[],
    datas=[
        (str(ROOT / "resources" / "bin"), "resources/bin"),
    ],
    hiddenimports=[
        "rdkit",
        "rdkit.Chem.AllChem",
        "meeko",
        "openbabel",
        "scipy.spatial",
        "gemmi",
        "uvicorn.lifespan",
        "uvicorn.protocols",
    ],
    hookspath=[],
    runtime_hooks=[],
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
