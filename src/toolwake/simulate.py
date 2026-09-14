"""Walk the path, lay the wake, and measure the tool against it.

The ordering inside the loop is the whole algorithm:

    for each row i
        1. place the tool at row i
        2. measure it against everything deposited before (i - lag)
        3. THEN deposit row i, if it is a print move

Checking before depositing is what stops the tool colliding with the bead it is
laying at that instant. `lag` extends that immunity a few rows back, because
material immediately behind the nozzle is still under the tool by construction.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .deposit import Deposit
from .tool import Needle
from .toolpath import PRINT, Toolpath

__all__ = ["Result", "simulate"]


@dataclass
class Result:
    """Per-row clearance, plus the wake that produced it."""

    path: Toolpath
    needle: Needle
    deposit: Deposit
    clearance: np.ndarray               # (N,) metres; inf where nothing in range
    culprit: np.ndarray                 # (N,) index of the closest bead, -1 if none
    source: np.ndarray                  # (N,) name of the closest tool part
    lag: int
    threshold: float
    meta: dict = field(default_factory=dict)

    @property
    def worst(self) -> float:
        finite = self.clearance[np.isfinite(self.clearance)]
        return float(finite.min()) if finite.size else float("inf")

    @property
    def n_hits(self) -> int:
        """Rows where the tool actually penetrates deposited material."""
        return int(np.sum(self.clearance < 0))

    @property
    def n_close(self) -> int:
        """Rows inside the warning band but not penetrating."""
        return int(np.sum((self.clearance >= 0) & (self.clearance < self.threshold)))

    @property
    def status(self) -> str:
        if self.n_hits:
            return "collision"
        if self.n_close:
            return "warn"
        return "ok"

    def report(self) -> dict:
        """A JSON-ready summary, shaped like a clearance report."""
        w = self.worst
        i = int(np.argmin(np.where(np.isfinite(self.clearance),
                                   self.clearance, np.inf)))
        return {
            "status": self.status,
            "rows": len(self.path),
            "deposited_segments": len(self.deposit),
            "min_clearance_mm": None if not np.isfinite(w) else round(w * 1e3, 4),
            "min_clearance_row": i if np.isfinite(w) else None,
            "min_clearance_source": str(self.source[i]) if np.isfinite(w) else None,
            "rows_penetrating": self.n_hits,
            "rows_below_threshold": self.n_close,
            "threshold_mm": round(self.threshold * 1e3, 4),
            "lag_rows": self.lag,
            **self.meta,
        }

    def __repr__(self):
        w = self.worst
        ws = "inf" if not np.isfinite(w) else f"{w*1e3:.3f} mm"
        return (f"Result({self.status}, {len(self.path)} rows, min={ws}, "
                f"{self.n_hits} penetrating)")


def simulate(path: Toolpath, needle: Needle | None = None, *, lag: int = 8,
             threshold: float = 5e-4, search: float = 0.02,
             progress=None) -> Result:
    """Run the wake simulation over `path`.

    Args:
        path: the toolpath to walk.
        needle: tool model; a default 34G-ish needle if omitted.
        lag: rows of immunity behind the nozzle. Material laid within this many
            rows is not treated as an obstacle.
        threshold: clearance below which a row counts as "close", metres.
        search: broad-phase radius, metres. Beads beyond this are not measured;
            raise it if the housing is large.
        progress: optional callable(i, n) for a progress bar.

    Returns:
        A `Result` carrying per-row clearance and the accumulated wake.
    """
    needle = needle or Needle()
    # The grid only exists to shrink the candidate set, so scale it to the
    # query, not the bead — see Deposit.__init__.
    dep = Deposit(bead_radius=needle.bead_radius,
                  cell=max(search / 4.0, 8.0 * needle.bead_radius))
    n = len(path)
    clearance = np.full(n, np.inf)
    culprit = np.full(n, -1, dtype=int)
    source = np.full(n, "", dtype=object)

    for i in range(n):
        rv = None if path.rotvec is None else path.rotvec[i]
        pose = needle.at(path.xyz[i], rv)

        best, who, which = np.inf, "", -1
        # Every section of the profile, plus the housing. Reporting WHICH part
        # was closest is the point: "hub" and "cannula" call for different
        # fixes, and a single "needle" label hides that.
        for name, shape in pose.named_shapes():
            d, j = dep.clearance(shape, frame=i, lag=lag, search=search)
            if d < best:
                best, who, which = d, name, j
        clearance[i] = best
        culprit[i] = which
        source[i] = who

        # Deposit AFTER measuring — see the module docstring.
        if i > 0 and path.kinds[i] == PRINT:
            dep.add(path.xyz[i - 1], path.xyz[i], frame=i)

        if progress is not None and (i % 200 == 0 or i == n - 1):
            progress(i + 1, n)

    return Result(path=path, needle=needle, deposit=dep, clearance=clearance,
                  culprit=culprit, source=source, lag=lag, threshold=threshold)
