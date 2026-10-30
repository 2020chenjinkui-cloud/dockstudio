"""Desktop sidecar lifecycle: parent watchdog + single-instance lock."""

from __future__ import annotations

import atexit
import os
import threading
import time
from pathlib import Path

_LOCK_PATH: Path | None = None
_WATCHDOG: threading.Thread | None = None


def _process_alive(pid: int) -> bool:
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes
        from ctypes import wintypes

        kernel32 = ctypes.windll.kernel32
        handle = kernel32.OpenProcess(0x1000, False, pid)
        if not handle:
            return False
        try:
            code = wintypes.DWORD()
            if not kernel32.GetExitCodeProcess(handle, ctypes.byref(code)):
                return True
            return code.value == 259  # STILL_ACTIVE
        finally:
            kernel32.CloseHandle(handle)
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True


def install_parent_watchdog(interval: float = 2.0) -> threading.Thread | None:
    """Exit the sidecar when its parent dies.

    PyInstaller onefile starts a bootloader parent plus an extracted child.  Tauri's
    CommandChild.kill() terminates the bootloader; without this watchdog the Python
    child keeps running and holds port 8000, so the next launch can never bind.
    """
    global _WATCHDOG
    if os.environ.get("DOCKSTUDIO_NO_PARENT_WATCHDOG") == "1":
        return None
    parent = os.getppid()
    if parent <= 1:
        return None
    if _WATCHDOG is not None and _WATCHDOG.is_alive():
        return _WATCHDOG

    def watch() -> None:
        while True:
            time.sleep(interval)
            if not _process_alive(parent):
                os._exit(0)

    _WATCHDOG = threading.Thread(target=watch, name="parent-watchdog", daemon=True)
    _WATCHDOG.start()
    return _WATCHDOG


def release_single_instance_lock() -> None:
    global _LOCK_PATH
    if _LOCK_PATH is None:
        return
    try:
        if _LOCK_PATH.read_text(encoding="utf-8").strip() == str(os.getpid()):
            _LOCK_PATH.unlink(missing_ok=True)
    except OSError:
        pass
    _LOCK_PATH = None


def acquire_single_instance_lock(data_dir: Path) -> bool:
    """Claim a PID file in the data dir; stale entries are reclaimed.

    A PID file avoids Windows byte-range locking, which is unreliable for the
    non-ASCII user paths this app targets and can leave a second launch spinning.
    """
    global _LOCK_PATH
    if os.environ.get("DOCKSTUDIO_NO_SINGLE_INSTANCE") == "1":
        return True
    path = Path(data_dir) / "backend.lock"
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return True
    for _ in range(2):
        try:
            with open(path, "x", encoding="utf-8") as handle:
                handle.write(str(os.getpid()))
            _LOCK_PATH = path
            atexit.register(release_single_instance_lock)
            return True
        except FileExistsError:
            pass
        except OSError:
            return True
        try:
            existing = int(path.read_text(encoding="utf-8").strip() or "0")
        except (OSError, ValueError):
            existing = 0
        if existing and _process_alive(existing):
            return False
        try:
            path.unlink(missing_ok=True)
        except OSError:
            return False
    return False
