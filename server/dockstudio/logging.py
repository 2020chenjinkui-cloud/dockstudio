"""Structured logging: console + rotating file under the data dir."""

from __future__ import annotations

import logging
from logging.handlers import RotatingFileHandler
from pathlib import Path

_FMT = "%(asctime)s %(levelname)s %(name)s: %(message)s"


def setup_logging(log_dir: Path, level: int = logging.INFO) -> None:
    log_dir.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger("dockstudio")
    root.setLevel(level)
    if root.handlers:
        return
    console = logging.StreamHandler()
    console.setFormatter(logging.Formatter(_FMT))
    file_handler = RotatingFileHandler(
        log_dir / "dockstudio.log", maxBytes=5_000_000, backupCount=3, encoding="utf-8"
    )
    file_handler.setFormatter(logging.Formatter(_FMT))
    root.addHandler(console)
    root.addHandler(file_handler)


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(f"dockstudio.{name}")
