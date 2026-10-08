"""PyInstaller runtime hook：为冻结的 openbabel 指定插件目录（BABEL_LIBDIR）。

openbabel 官方轮子的格式插件布局（实测）：
- macOS/Linux: openbabel/lib/openbabel/<ver>/*.so
- Windows:     openbabel/bin/*.obf（3.1.x openbabel-wheel 在 lib/）
冻结后收集目录在 sys._MEIPASS 下，按含格式插件文件的目录探测。
"""

import os
import sys
from pathlib import Path

_meipass = Path(getattr(sys, "_MEIPASS", ""))
_PLUGIN_SUFFIXES = {".so", ".dll", ".obf"}


def _looks_like_plugin_dir(d: Path) -> bool:
    try:
        return any(
            f.is_file() and f.suffix.lower() in _PLUGIN_SUFFIXES and "format" in f.name.lower()
            for f in d.iterdir()
        )
    except OSError:
        return False


def _find_plugin_dir() -> Path | None:
    candidates = [
        _meipass / "openbabel" / "lib" / "openbabel",  # mac/linux: <ver>/*.so 的上一级
        _meipass / "_internal" / "openbabel" / "lib" / "openbabel",
        _meipass / "openbabel" / "bin",
        _meipass / "_internal" / "openbabel" / "bin",
        _meipass / "openbabel" / "lib",
        _meipass / "_internal" / "openbabel" / "lib",
    ]
    expanded: list[Path] = []
    for c in candidates:
        if c.is_dir():
            expanded.append(c)
            try:
                expanded.extend([d for d in c.iterdir() if d.is_dir()])
            except OSError:
                continue
    for d in expanded:
        if _looks_like_plugin_dir(d):
            return d
    # 兜底：全局搜含 format 插件文件的目录
    for root in (_meipass, _meipass / "_internal"):
        if not root.is_dir():
            continue
        for f in root.rglob("*format*"):
            if f.is_file() and f.suffix.lower() in _PLUGIN_SUFFIXES:
                return f.parent
    return None


_plugin_dir = _find_plugin_dir()
if _plugin_dir is not None:
    os.environ.setdefault("BABEL_LIBDIR", str(_plugin_dir))

# 数据文件（力场参数等）位置：openbabel/share/openbabel/<ver>
for _cand in (
    _meipass / "openbabel" / "share" / "openbabel",
    _meipass / "_internal" / "openbabel" / "share" / "openbabel",
):
    if _cand.is_dir():
        for _ver in _cand.iterdir():
            if _ver.is_dir():
                os.environ.setdefault("BABEL_DATADIR", str(_ver))
                break
        break
