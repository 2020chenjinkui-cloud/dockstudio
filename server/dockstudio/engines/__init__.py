"""Engine registry: discover adapters and expose them to the API layer."""

from __future__ import annotations

from .base import DockEngine

__all__ = ["available_engines", "get_engine"]


def available_engines() -> list[DockEngine]:
    """Instantiate all built-in engine adapters. Unavailable engines are still
    returned (manifest says why) so the UI can grey them out."""
    from .smina import SminaEngine
    from .vina import VinaEngine

    return [VinaEngine(), SminaEngine()]


def get_engine(name: str) -> DockEngine | None:
    for engine in available_engines():
        if engine.manifest.name == name:
            return engine
    return None
