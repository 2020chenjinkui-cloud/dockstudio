"""PyInstaller runtime hook：为冻结的 openbabel 指定插件目录（BABEL_LIBDIR）。"""

import os
import sys
from pathlib import Path

_meipass = Path(getattr(sys, "_MEIPASS", ""))

def _set_if_present(*parts: str, env: str) -> None:
    for root in (_meipass, _meipass / "openbabel", _meipass / "_internal"):
        candidate = root.joinpath(*parts)
        if candidate.is_dir():
            os.environ.setdefault(env, str(candidate))
            return


# openbabel-wheel 的格式插件通常随 collect_all 落在 <meipass>/openbabel/lib
_set_if_present("lib", env="BABEL_LIBDIR")
_set_if_present("openbabel", "lib", env="BABEL_LIBDIR")
if "BABEL_LIBDIR" not in os.environ:
    candidate = _meipass / "openbabel"
    if candidate.is_dir():
        os.environ.setdefault("BABEL_LIBDIR", str(candidate))
