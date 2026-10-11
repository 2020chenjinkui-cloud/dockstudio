"""EventLog: write/read/filter, blank-line tolerance, multi-emit."""

from __future__ import annotations

import json
from pathlib import Path

from dockstudio.core.events import EventLog


def test_emit_and_read(tmp_path: Path) -> None:
    log = EventLog(tmp_path)
    r1 = log.emit("job.created", "job", "j1", {"engine": "vina", "n": 1})
    log.emit("job.completed", "job", "j1", {"top": -9.2})
    log.emit("receptor.prepared", "receptor", "r1", {"chains": 2}, actor="system")
    assert r1["subject"] == {"kind": "job", "id": "j1"}
    assert r1["actor"] == "user"
    events = list(log.read())
    assert len(events) == 3
    assert events[0]["action"] == "job.created"
    assert events[1]["payload"]["top"] == -9.2
    assert events[2]["actor"] == "system"
    # jsonl is one object per line, ISO8601 ts
    lines = (tmp_path / "events.jsonl").read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 3
    assert all(json.loads(line)["ts"] for line in lines)


def test_read_filters(tmp_path: Path) -> None:
    log = EventLog(tmp_path)
    log.emit("a1", "job", "j1", {})
    log.emit("a2", "job", "j2", {})
    log.emit("a3", "ligand", "l1", {})
    assert [e["action"] for e in log.read(filter_kind="job")] == ["a1", "a2"]
    assert [e["action"] for e in log.read(filter_id="j2")] == ["a2"]
    since = next(iter(log.read()))["ts"]
    # since returns events strictly newer than the given ts
    newer = [e for e in log.read(since=since)]
    assert all(e["action"] != "a1" for e in newer)


def test_blank_and_garbage_lines_tolerated(tmp_path: Path) -> None:
    log = EventLog(tmp_path)
    log.emit("a1", "job", "j1", {})
    # simulate a crash-torn tail: partial line + blank lines
    with (tmp_path / "events.jsonl").open("a", encoding="utf-8") as fh:
        fh.write("\n{'partial json\n\n")
    log.emit("a2", "job", "j2", {})
    actions = [e["action"] for e in log.read()]
    assert actions == ["a1", "a2"]


def test_read_missing_file(tmp_path: Path) -> None:
    log = EventLog(tmp_path / "nonexistent")
    assert list(log.read()) == []
