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
             bead_radius: float | None = None, bead_drop: float = 0.0,
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
        bead_radius: radius of the laid bead, metres. Defaults to the needle's
            bore radius, which assumes the bead keeps the round cross-section
            it had inside the needle.

            It does not. A bead laid at layer height h is squashed to h and
            spreads sideways, and modelling it round makes it TALLER than the
            layer it lives in — so it pokes up through where the nozzle will
            sit on the next pass and reads as a collision on rows that printed
            perfectly. On a 0.036 mm-layer slab with a 0.09 mm bore that was
            498 rows reported, 391 of them blamed on material one layer below.
            Passing `min(bore_radius, h / 2)` leaves 41, and those are genuine
            same-layer retraces.

            A bead cannot be taller than the layer it was laid in. Left as an
            explicit argument rather than inferred, because only the caller
            knows the process.
        bead_drop: how far BELOW the nozzle face the bead sits, metres,
            measured along the tool axis. Zero puts the bead's centre on the
            toolpath, which is where the nozzle tip is — so every bead ends up
            half-embedded in the plane the nozzle face travels in, and any
            same-layer neighbour passing under the wall reports a collision of
            exactly one bead radius. Those are tangent contacts dressed up as
            penetrations: on the reference slab, 41 rows, every one of them at
            exactly -bead_radius.

            Material leaving the bore fills the gap between the previous
            layer's top and the nozzle face, so it occupies [z - h, z] and its
            centre is h/2 down. Passing `h / 2` alongside `bead_radius=h / 2`
            makes the bead exactly fill its layer, tangent to the face that
            laid it. That takes the same slab to zero.
        progress: optional callable(i, n) for a progress bar.

    Returns:
        A `Result` carrying per-row clearance and the accumulated wake.
    """
    needle = needle or Needle()
    r_bead = needle.bead_radius if bead_radius is None else float(bead_radius)
    if r_bead < 0:
        raise ValueError("bead_radius must be >= 0")
    # The grid only exists to shrink the candidate set, so scale it to the
    # query, not the bead — see Deposit.__init__.
    dep = Deposit(bead_radius=r_bead,
                  cell=max(search / 4.0, 8.0 * r_bead))
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
            a, b = path.xyz[i - 1], path.xyz[i]
            if bead_drop:
                # Along the TOOL axis, not -Z: on a non-planar move the
                # material still lands under the nozzle face, wherever that
                # face is pointing.
                shift = pose.axis * float(bead_drop)
                a, b = a - shift, b - shift
            dep.add(a, b, frame=i)

        if progress is not None and (i % 200 == 0 or i == n - 1):
            progress(i + 1, n)

    return Result(path=path, needle=needle, deposit=dep, clearance=clearance,
                  culprit=culprit, source=source, lag=lag, threshold=threshold)
