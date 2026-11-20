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

Two grids, because the queries come in two sizes. A needle query wants cells
about the size of the needle; a hub or housing query wants cells about the size
of the hub. One grid sized for either makes the other pay — see `_candidates`.
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["Deposit"]

# Starting capacity. Small enough not to matter, big enough that short
# toolpaths never reallocate.
_INITIAL_CAPACITY = 1024

# Most fine cells a query may inspect before it moves to the coarse grid.
# Above this a query is hub- or housing-sized, and walking it at needle
# resolution costs more than the extra candidates the coarse grid hands back.
_FINE_PROBE_LIMIT = 512


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
        coarse_cell: cell size of a second, coarser grid that queries too big
            for the fine one fall back to, metres. Defaults to `cell`, which
            means one grid. `simulate` sizes `cell` to the needle tip and this
            to the widest part of the tool.
    """

    def __init__(self, bead_radius: float, cell: float | None = None,
                 coarse_cell: float | None = None):
        if bead_radius < 0:
            raise ValueError("bead_radius must be >= 0")
        self.bead_radius = float(bead_radius)
        self.cell = float(cell) if cell else max(8.0 * self.bead_radius, 5e-3)
        self.coarse_cell = float(coarse_cell) if coarse_cell else self.cell
        if self.coarse_cell < self.cell:
            raise ValueError("coarse_cell must be >= cell")
        # Stored as capacity-doubling arrays, NOT Python lists.
        #
        # They were lists, and `clearance` did `np.asarray(self._a)[cand]` to
        # index them. That rebuilds an array of the ENTIRE wake on every query,
        # so each row cost O(total deposit) no matter how few candidates the
        # spatial hash returned — making the whole sweep quadratic and hiding
        # it behind a broad phase that was doing its job perfectly. Measured
        # with the candidate count pinned at 2, a query went 673 us at 500
        # segments to 12.5 ms at 16 000: pure deposit-size cost, zero of it
        # real work.
        cap = _INITIAL_CAPACITY
        self._n = 0
        self._A = np.empty((cap, 3), dtype=float)   # segment starts
        self._B = np.empty((cap, 3), dtype=float)   # segment ends
        self._T = np.empty(cap, dtype=np.int64)     # frame each was laid at
        self._grid: dict[tuple, list[int]] = {}
        self._coarse: dict[tuple, list[int]] | None = (
            {} if self.coarse_cell > self.cell else None)
        # Bounding box of every measured point laid so far, padded by the bead
        # radius. A query that misses it cannot hit anything, which on a planar
        # print is the hub, collar and housing on every row: they sit above
        # everything printed. Python floats, because this is checked per query.
        self._lo = [float("inf")] * 3
        self._hi = [float("-inf")] * 3
        # Cells inspected by the most recent _candidates call. Exposed so the
        # cost invariant below can be tested without timing anything.
        self.last_probe = 0

    # ---------------------------------------------------------------- build
    def __len__(self) -> int:
        return self._n

    # Views over the live prefix. Kept under the old names so callers that
    # reach in (viewer.py reads _t) keep working; they are now O(1) slices
    # rather than list-to-array conversions.
    @property
    def _a(self) -> np.ndarray:
        return self._A[:self._n]

    @property
    def _b(self) -> np.ndarray:
        return self._B[:self._n]

    @property
    def _t(self) -> np.ndarray:
        return self._T[:self._n]

    def _ensure(self, n: int) -> None:
        """Grow to hold at least `n` segments, doubling so `add` stays O(1)
        amortised."""
        cap = self._T.shape[0]
        if n <= cap:
            return
        cap = max(cap * 2, n)
        for name, shape in (("_A", (cap, 3)), ("_B", (cap, 3)), ("_T", (cap,))):
            old_arr = getattr(self, name)
            new_arr = np.empty(shape, dtype=old_arr.dtype)
            new_arr[:self._n] = old_arr[:self._n]
            setattr(self, name, new_arr)

    def _cells_for(self, a: np.ndarray, b: np.ndarray, cell: float) -> set:
        """Grid cells holding the points `clearance` measures on this bead.

        Those are the two ends and the midpoint — the narrow phase looks at
        nothing else — so a bead is filed under the cells of those three
        points, each padded by the bead radius.

        It used to be filed under every cell of the segment's whole bounding
        box. On a long diagonal infill line that is quadratic in its length: a
        200 mm diagonal at 2.75 mm cells is 5 000 cells, at needle-sized cells
        it is 150 000, for a bead whose measured points occupy at most 24.
        """
        pts = np.array([a, b, 0.5 * (a + b)])
        lo = np.floor((pts - self.bead_radius) / cell).astype(np.int64).tolist()
        hi = np.floor((pts + self.bead_radius) / cell).astype(np.int64).tolist()
        out = set()
        for (l0, l1, l2), (h0, h1, h2) in zip(lo, hi):
            for i in range(l0, h0 + 1):
                for j in range(l1, h1 + 1):
                    for k in range(l2, h2 + 1):
                        out.add((i, j, k))
        return out

    def add(self, a, b, frame: int) -> None:
        """Record one deposited segment, laid at `frame`."""
        self._ensure(self._n + 1)
        idx = self._n
        self._A[idx] = a
        self._B[idx] = b
        self._T[idx] = frame
        self._n = idx + 1
        A, B = self._A[idx], self._B[idx]
        r = self.bead_radius
        for ax, (u, v) in enumerate(zip(A.tolist(), B.tolist())):
            self._lo[ax] = min(self._lo[ax], u - r, v - r)
            self._hi[ax] = max(self._hi[ax], u + r, v + r)
        for c in self._cells_for(A, B, self.cell):
            self._grid.setdefault(c, []).append(idx)
        if self._coarse is not None:
            for c in self._cells_for(A, B, self.coarse_cell):
                self._coarse.setdefault(c, []).append(idx)

    # ---------------------------------------------------------------- query
    def _candidates(self, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
        """Segment indices whose cells overlap the box [lo, hi].

        Two ways to answer this, and which is cheaper depends entirely on the
        query. Enumerating the box's cells costs one dict lookup per cell;
        scanning the occupied cells costs one box test per cell that actually
        holds material. Space is mostly empty, so for a large query the second
        is enormously cheaper — and tool queries ARE large, because the tool's
        bounding box includes the 130 mm luer housing while the grid is sized
        for a 45 um bead.

        Measured on small_slab.gcode before this branch existed: 43.8 million
        dict lookups over 8 715 queries, about 5 000 cells walked per query to
        find a few hundred occupied ones, and 86% of total runtime spent right
        here. Picking the cheaper enumeration is not an optimisation of the
        algorithm; it is declining to walk empty space.

        Both branches return the same set, so this is a pure cost choice.

        That choice stops being enough once the grid is fine enough to suit
        the needle. The occupied-cell scan is then O(material laid so far) on
        every hub and housing query, every row — a quadratic sweep again,
        measured on a 10 mm part as 336 -> 536 cells scanned per query over
        12 000 rows. So a query that would inspect more than
        `_FINE_PROBE_LIMIT` fine cells goes to the coarse grid instead, where
        a hub-sized box is a few dozen cells.

        The coarse grid is used when the box is at least one coarse cell wide
        on every axis — rounding it out to coarse cells then at most triples
        each side — or when the fine walk would be long anyway. Keying it only
        on the walk length left small parts on the fine grid, scanning up to
        `_FINE_PROBE_LIMIT` occupied cells for every hub query.
        """
        lo, hi = lo.tolist(), hi.tolist()
        grid, cell = self._grid, self.cell
        lo_i, hi_i, n_cells = self._box(lo, hi, cell)
        if self._coarse is not None and (
                min(h - l for l, h in zip(lo, hi)) >= self.coarse_cell
                or min(n_cells, len(grid)) > _FINE_PROBE_LIMIT):
            grid, cell = self._coarse, self.coarse_cell
            lo_i, hi_i, n_cells = self._box(lo, hi, cell)
        out: set[int] = set()

        self.last_probe = min(n_cells, len(grid))
        if n_cells <= len(grid):
            for i in range(lo_i[0], hi_i[0] + 1):
                for j in range(lo_i[1], hi_i[1] + 1):
                    for k in range(lo_i[2], hi_i[2] + 1):
                        hit = grid.get((i, j, k))
                        if hit:
                            out.update(hit)
        else:
            lo0, lo1, lo2 = lo_i
            hi0, hi1, hi2 = hi_i
            for (i, j, k), hit in grid.items():
                if lo0 <= i <= hi0 and lo1 <= j <= hi1 and lo2 <= k <= hi2:
                    out.update(hit)

        return np.fromiter(out, dtype=int, count=len(out))

    @staticmethod
    def _box(lo: list, hi: list, cell: float):
        """Cell index range of the box [lo, hi], and how many cells it spans."""
        lo_i = [math.floor(x / cell) for x in lo]
        hi_i = [math.floor(x / cell) for x in hi]
        # Python ints: the product of three cell counts overflows int64 for a
        # degenerate query, and an overflowed negative would silently pick the
        # wrong branch.
        n_cells = ((hi_i[0] - lo_i[0] + 1) * (hi_i[1] - lo_i[1] + 1)
                   * (hi_i[2] - lo_i[2] + 1))
        return lo_i, hi_i, n_cells

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
            search: broad-phase radius, metres. Every bead whose measured
                points come within `search` of the shape is measured, so a
                result <= `search` is exact. A result above it is only
                whatever the grid happened to hand back — the nearest bead may
                not have been looked at — and is not a distance to rely on.
                `simulate` widens the search until it gets an exact answer, or
                reports `inf`.

        Returns:
            (clearance, segment_index). `inf` and -1 when nothing is in range.
        """
        if self._n == 0:
            return float("inf"), -1

        lo, hi = shape.bounds(pad=search)
        if any(l > H or h < L for l, h, L, H
               in zip(lo.tolist(), hi.tolist(), self._lo, self._hi)):
            return float("inf"), -1
        cand = self._candidates(lo, hi)
        if cand.size == 0:
            return float("inf"), -1

        # Fancy-indexing the live arrays costs O(candidates). The old
        # np.asarray(list)[cand] cost O(total deposit) on every single query.
        cand = cand[self._T[cand] < frame - lag]   # earlier material only
        if cand.size == 0:
            return float("inf"), -1

        A = self._A[cand]
        B = self._B[cand]
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
        if self._n == 0:
            return np.zeros((0, 2, 3))
        A = self._A[:self._n]
        B = self._B[:self._n]
        if upto is not None:
            keep = self._T[:self._n] <= upto
            A, B = A[keep], B[keep]
        return np.stack([A, B], axis=1)

    def bounds(self):
        if self._n == 0:
            return None
        pts = np.vstack([self._A[:self._n], self._B[:self._n]])
        return pts.min(axis=0) - self.bead_radius, pts.max(axis=0) + self.bead_radius

    def __repr__(self):
        return (f"Deposit({len(self)} segments, bead_r={self.bead_radius:.4g}, "
                f"cell={self.cell:.4g}, {len(self._grid)} cells)")
