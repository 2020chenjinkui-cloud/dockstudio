"""JSONL event stream: append-only audit log per project directory."""

from __future__ import annotations

import json
import threading
from collections.abc import Iterator
from datetime import datetime, timezone
from pathlib import Path

_LOCKS: dict[str, threading.Lock] = {}
_LOCKS_GUARD = threading.Lock()


def _lock_for(path: Path) -> threading.Lock:
    key = str(path)
    with _LOCKS_GUARD:
        return _LOCKS.setdefault(key, threading.Lock())


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds")


class EventLog:
    """Append-only JSONL event log scoped to one project directory."""

    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir)
        self.path = self.project_dir / "events.jsonl"

    def emit(
        self,
        action: str,
        subject_kind: str,
        subject_id: str,
        payload: dict,
        actor: str = "user",
    ) -> dict:
        """Append one event line. Thread- and process-safe via per-path lock;
        the OS-level O_APPEND keeps same-host concurrent writers from tearing
        lines for the small payloads we use."""
        record = {
            "ts": _utc_now(),
            "actor": actor,
            "action": action,
            "subject": {"kind": subject_kind, "id": subject_id},
            "payload": payload,
        }
        line = json.dumps(record, ensure_ascii=False, separators=(",", ":"))
        lock = _lock_for(self.path)
        with lock:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            with self.path.open("a", encoding="utf-8") as fh:
                fh.write(line + "\n")
        return record

    def read(
        self,
        filter_kind: str | None = None,
        filter_id: str | None = None,
        since: str | None = None,
    ) -> Iterator[dict]:
        """Yield events, tolerating blank/garbled lines (crash-torn tail)."""
        if not self.path.is_file():
            return
        lock = _lock_for(self.path)
        with lock:
            text = self.path.read_text(encoding="utf-8", errors="replace")
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            subject = event.get("subject") or {}
            if filter_kind is not None and subject.get("kind") != filter_kind:
                continue
            if filter_id is not None and subject.get("id") != filter_id:
                continue
            if since is not None and str(event.get("ts", "")) <= since:
                continue
            yield event


def sha256_file(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def env_fingerprint():
    """Snapshot the software environment for job.json (reproducibility)."""
    import platform
    from importlib.metadata import PackageNotFoundError, version

    from ..models import EnvFingerprint

    def ver(name: str) -> str:
        try:
            return version(name)
        except PackageNotFoundError:
            return ""

    return EnvFingerprint(
        python=platform.python_version(),
        rdkit=ver("rdkit"),
        meeko=ver("meeko"),
        openbabel=ver("openbabel-wheel") or ver("openbabel"),
        platform=f"{platform.system()} {platform.release()} {platform.machine()}",
    )
