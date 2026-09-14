"""The wake: material already laid down, and how to ask what is near it.

Two-phase by design.

  Broad phase  — a uniform spatial hash (a voxel grid holding segment indices).
                 O(1) lookup, answers "which beads could possibly be near here".
  Narrow phase — exact capsule distance to those candidate beads only.

A plain occupancy grid would give a binary answer quantised to the voxel size;
keeping the segments and measuring against them means the clearance number in
the report is a real distance, not a voxel count. The grid exists purely to
avoid testing every bead against every tool pose.

Time ordering is the other reason this class exists. Material can only obstruct
the tool if it was deposited EARLIER, and not so recently that it is still
under the nozzle. `add` appends in path order and `query` takes the frame index,
so the trailing-window rule lives in one place instead of at every call site.
"""
from __future__ import annotations

import numpy as np

__all__ = ["Deposit"]


class Deposit:
    """Accumulates deposited bead segments and answers clearance queries.

    Args:
        bead_radius: radius of the laid bead, metres. For a well-tuned print
            this is about half the needle inner diameter.
        cell: spatial-hash cell size, metres. This wants to be on the order of
            the QUERY radius, not the bead. Sized to the bead (45 um) against a
            20 mm search, a single lookup walks 41^3 = 69k cells and the sweep
            crawls; at 5 mm it walks 9^3 = 729. Callers that know their search
            radius should pass `search / 4`.
    """

    def __init__(self, bead_radius: float, cell: float | None = None):
        if bead_radius < 0:
            raise ValueError("bead_radius must be >= 0")
        self.bead_radius = float(bead_radius)
        self.cell = float(cell) if cell else max(8.0 * self.bead_radius, 5e-3)
        self._a: list[np.ndarray] = []       # segment starts
        self._b: list[np.ndarray] = []       # segment ends
        self._t: list[int] = []              # frame index each was laid at
        self._grid: dict[tuple, list[int]] = {}

    # ---------------------------------------------------------------- build
    def __len__(self) -> int:
        return len(self._a)

    def _cells_for(self, a: np.ndarray, b: np.ndarray):
        """Every grid cell the segment's padded AABB touches."""
        lo = np.minimum(a, b) - self.bead_radius
        hi = np.maximum(a, b) + self.bead_radius
        lo_i = np.floor(lo / self.cell).astype(int)
        hi_i = np.floor(hi / self.cell).astype(int)
        for i in range(lo_i[0], hi_i[0] + 1):
            for j in range(lo_i[1], hi_i[1] + 1):
                for k in range(lo_i[2], hi_i[2] + 1):
                    yield (i, j, k)

    def add(self, a, b, frame: int) -> None:
        """Record one deposited segment, laid at `frame`."""
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        idx = len(self._a)
        self._a.append(a)
        self._b.append(b)
        self._t.append(int(frame))
        for c in self._cells_for(a, b):
            self._grid.setdefault(c, []).append(idx)

    # ---------------------------------------------------------------- query
    def _candidates(self, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        lo_i = np.floor(lo / self.cell).astype(int)
        hi_i = np.floor(hi / self.cell).astype(int)
        out: set[int] = set()
        for i in range(lo_i[0], hi_i[0] + 1):
            for j in range(lo_i[1], hi_i[1] + 1):
                for k in range(lo_i[2], hi_i[2] + 1):
                    hit = self._grid.get((i, j, k))
                    if hit:
                        out.update(hit)
        return np.fromiter(out, dtype=int, count=len(out))

    def clearance(self, shape, frame: int, lag: int = 0,
                  search: float = 0.02) -> tuple[float, int]:
        """Smallest gap between `shape` and material laid before `frame - lag`.

        Args:
            shape: anything with `.distance(points)` and `.bounds(pad)` —
                a Capsule or a Box from `toolwake.geometry`.
            frame: the current frame index.
            lag: how many frames back to start counting material as an
                obstacle. Without this the tool always collides with the bead
                it is extruding right now.
            search: broad-phase radius, metres. Beads further away than this
                are not measured at all; the returned distance is clipped to it.

        Returns:
            (clearance, segment_index). `inf` and -1 when nothing is in range.
        """
        if not self._a:
            return float("inf"), -1

        lo, hi = shape.bounds(pad=search)
        cand = self._candidates(lo, hi)
        if cand.size == 0:
            return float("inf"), -1

        t = np.asarray(self._t)[cand]
        cand = cand[t < frame - lag]              # earlier material only
        if cand.size == 0:
            return float("inf"), -1

        A = np.asarray(self._a)[cand]
        B = np.asarray(self._b)[cand]
        # Sample each candidate bead along its axis; with beads this short
        # relative to the tool, endpoints plus midpoint bound the true
        # segment-to-segment distance closely and stay fully vectorised.
        pts = np.vstack([A, B, 0.5 * (A + B)])
        d = shape.distance(pts) - self.bead_radius
        n = len(cand)
        d = np.min(d.reshape(3, n), axis=0)
        j = int(np.argmin(d))
        return float(d[j]), int(cand[j])

    # ---------------------------------------------------------------- views
    def segments(self, upto: int | None = None) -> np.ndarray:
        """(N, 2, 3) array of laid segments, for drawing. `upto` filters by frame."""
        if not self._a:
            return np.zeros((0, 2, 3))
        A = np.asarray(self._a)
        B = np.asarray(self._b)
        if upto is not None:
            keep = np.asarray(self._t) <= upto
            A, B = A[keep], B[keep]
        return np.stack([A, B], axis=1)

    def bounds(self):
        if not self._a:
            return None
        pts = np.vstack([np.asarray(self._a), np.asarray(self._b)])
        return pts.min(axis=0) - self.bead_radius, pts.max(axis=0) + self.bead_radius

    def __repr__(self):
        return (f"Deposit({len(self)} segments, bead_r={self.bead_radius:.4g}, "
                f"cell={self.cell:.4g}, {len(self._grid)} cells)")
