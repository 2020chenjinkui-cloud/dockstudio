"""PyInstaller runtime hook: make RDKit data readable on non-ASCII Windows paths.

RDKit's C++ code opens Data/BaseFeatures.fdef with a narrow-char path.  When the
onefile temp directory contains a non-ASCII user name, PoseBusters and feature
factories fail at import.  Copy the small Data directory to an ASCII location and
set RDBASE before the application imports rdkit.
"""

import os
import shutil
import sys
from pathlib import Path


def _ensure() -> None:
    if os.environ.get("RDBASE"):
        return
    meipass = Path(getattr(sys, "_MEIPASS", ""))
    src = meipass / "rdkit" / "Data"
    if not src.is_dir() or all(ord(ch) < 128 for ch in str(src)):
        return
    for base in (os.environ.get("PUBLIC"), os.environ.get("PROGRAMDATA")):
        if not base:
            continue
        target = Path(base) / "DockStudio" / "rdkit-data"
        if not all(ord(ch) < 128 for ch in str(target)):
            continue
        try:
            target.mkdir(parents=True, exist_ok=True)
            probe = target / "Data" / "BaseFeatures.fdef"
            if not probe.is_file():
                shutil.copytree(src, target / "Data", dirs_exist_ok=True)
            if probe.is_file():
                os.environ["RDBASE"] = str(target)
                return
        except OSError:
            continue


_ensure()
