"""Degradation ledger: every silent fallback becomes a queryable record.

The review's recurring criticism was that DockStudio degrades quietly -- PoseBusters
falls back to 3 rules, gnina rescore fails into an empty column, a co-crystal SDF has
no bond orders -- and the report never says so.  Rather than chase 165 individual
``except`` blocks, this gives them one place to report to, so a run can answer
"what did I not actually measure?" in one call.
"""

from __future__ import annotations

import threading
from collections import Counter
from dataclasses import dataclass, field
from typing import Any


@dataclass(frozen=True)
class Degradation:
    """One thing the pipeline could not do as advertised."""

    category: str
    detail: str
    #: Where it happened, as "module.function" -- enough to find it, not a stack.
    origin: str = ""
    #: True when the affected number feeds a scientific conclusion (score, RMSD,
    #: funnel count).  Those must be surfaced to the user; cosmetic ones need not be.
    affects_science: bool = True
    extra: dict[str, Any] = field(default_factory=dict)

    def marker(self) -> str:
        """Stable, greppable token for CSV/JSON/audit consumers."""
        return f"DEGRADED {self.category}: {self.detail}"


class DegradationLedger:
    """Thread-safe, process-wide record of degradations for one run."""

    def __init__(self) -> None:
        self._items: list[Degradation] = []
        self._lock = threading.Lock()

    def record(
        self,
        category: str,
        detail: str,
        *,
        origin: str = "",
        affects_science: bool = True,
        **extra: Any,
    ) -> Degradation:
        item = Degradation(
            category=category,
            detail=detail,
            origin=origin,
            affects_science=affects_science,
            extra=dict(extra),
        )
        with self._lock:
            self._items.append(item)
        return item

    def items(self) -> list[Degradation]:
        with self._lock:
            return list(self._items)

    def science_items(self) -> list[Degradation]:
        return [i for i in self.items() if i.affects_science]

    def counts(self) -> dict[str, int]:
        return dict(Counter(i.category for i in self.items()))

    def is_clean(self) -> bool:
        """True only when nothing degraded -- the precondition for claiming a full run."""
        return not self._items

    def summary(self) -> str:
        if not self._items:
            return "no degradation recorded"
        parts = [f"{cat}×{n}" for cat, n in sorted(self.counts().items())]
        return "degraded: " + ", ".join(parts)

    def as_rows(self) -> list[dict[str, Any]]:
        return [
            {
                "category": i.category,
                "detail": i.detail,
                "origin": i.origin,
                "affects_science": i.affects_science,
                "marker": i.marker(),
                **i.extra,
            }
            for i in self.items()
        ]

    def clear(self) -> None:
        with self._lock:
            self._items.clear()


#: Process-wide default ledger.  A run that wants isolation can construct its own.
LEDGER = DegradationLedger()
