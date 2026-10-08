"""PyInstaller runtime hook：为冻结的 openbabel 指定插件目录（BABEL_LIBDIR）。

openbabel 3.2.x 官方轮子的格式插件位于包内 bin/（*.obf），
旧 openbabel-wheel 在 lib/。冻结后收集目录在 sys._MEIPASS 下，
按含 *.obf 的目录搜索最稳。
"""

import os
import sys
from pathlib import Path

_meipass = Path(getattr(sys, "_MEIPASS", ""))


def _find_plugin_dir() -> Path | None:
    candidates = [
        _meipass / "openbabel" / "bin",
        _meipass / "openbabel" / "lib",
        _meipass / "openbabel",
        _meipass / "_internal" / "openbabel" / "bin",
        _meipass / "_internal" / "openbabel" / "lib",
    ]
    for c in candidates:
        if c.is_dir():
            try:
                if any(f.suffix.lower() in (".obf", ".so", ".dll") and "format" in f.name.lower() for f in c.iterdir() if f.is_file()):
                    return c
            except OSError:
                continue
    # 兜底：一层搜含 .obf 的目录
    for root in (_meipass, _meipass / "_internal"):
        if not root.is_dir():
            continue
        for d in root.rglob("*.obf"):
            return d.parent
    return None


_plugin_dir = _find_plugin_dir()
if _plugin_dir is not None:
    os.environ.setdefault("BABEL_LIBDIR", str(_plugin_dir))
