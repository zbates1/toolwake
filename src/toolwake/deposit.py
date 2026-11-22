"""The wake: material already laid down, and how to ask what is near it.

Two-phase by design.

  Broad phase  — a uniform spatial hash (a voxel grid holding segment indices).
                 O(1) lookup, answers "which beads could possibly be near here".
  Narrow phase — exact distance from the tool to the closest point along each
                 candidate bead, found by branch and bound (`_min_along`).

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

# Beads longer than this many fine cells are filed in the coarse grid only.
# See `Deposit.add`.
_LONG_BEAD_CELLS = 8

# How close to the true minimum along a bead the narrow phase proves it is,
# metres. The reported clearance is never below the truth and at most this far
# above it. A nanometre is a thousandth of anything the bead model resolves.
_NARROW_TOL = 1e-9

# Ceilings on the narrow phase's search, so a pathological shape cannot run
# away with memory or time. Neither is reached by the shapes in geometry.py on
# any toolpath tested; if one is, the answer is still a real measured
# distance, just not one proven to be within _NARROW_TOL of the minimum.
_MAX_ROUNDS = 200
_MAX_INTERVALS = 1 << 16


def _min_along(shape, A: np.ndarray, B: np.ndarray, cutoff: float,
               tol: float = _NARROW_TOL) -> tuple[float, int]:
    """Smallest distance from `shape` to any point on the segments A[k]-B[k].

    Branch and bound along each segment. A stretch of segment is split until
    it is proven unable to beat the best distance measured so far by more
    than `tol`, using two lower bounds on what it could contain:

      Lipschitz  a distance field changes no faster than the point moves, so
                 a stretch of length l with ends at f0 and f1 reaches no lower
                 than (f0 + f1 - l) / 2.
      Convex     a convex function never dips below its tangent lines. Shapes
                 that provide `_probe` give a convex lower bound on their
                 distance and its gradient, which pins a smooth minimum in a
                 split or two and a flat one immediately — where the
                 Lipschitz bound alone would split a bead lying flat under
                 the nozzle face down to the nanometre.

    Each segment is measured first at its ends, then — if those cannot rule it
    out — wherever the shape's `_hints` say its closest point is likely, and
    after that at the best guesses `_bound` makes from the ends of each
    stretch still alive. On the toolpaths tested most queries are settled by
    the hints alone.

    Anything that cannot come within `cutoff` is never refined. Returns the
    smallest distance measured and the segment it was measured on: a real
    distance to a real point, never below the true minimum, and within `tol`
    of it wherever that minimum is within `cutoff`.

    This replaces measuring each bead at its two ends and midpoint, which was
    exact for beads short against the tool and blind to the middle of long
    ones: a needle pushed into a 10 mm bead 2.5 mm from one end measured the
    bead as 2.5 mm away.
    """
    probe = getattr(shape, "_probe", None)
    if probe is None:
        def probe(P):
            # No convex bound: one that is never informative.
            f = shape.distance(P)
            return f, np.full_like(f, -np.inf), np.zeros_like(P)

    n = len(A)
    D = B - A
    f, h, g = probe(np.concatenate([A, B]))
    s = np.einsum("ij,ij->i", g, np.concatenate([D, D]))
    ends = np.minimum(f[:n], f[n:])
    j = int(np.argmin(ends))
    best, who = float(ends[j]), j

    # One row per live stretch of segment. Columns: the stretch's ends in
    # segment parameter (0 at A, 1 at B), the distance at each end, the convex
    # bound and its slope along the segment at each end, then the segment's
    # start, direction, length and index. One array so that keeping and
    # splitting stretches is one operation, not one per quantity.
    st = np.empty((n, 16))
    st[:, 0], st[:, 1] = 0.0, 1.0
    st[:, 2], st[:, 3] = f[:n], f[n:]
    st[:, 4], st[:, 5] = h[:n], h[n:]
    st[:, 6], st[:, 7] = s[:n], s[n:]
    st[:, 8:11], st[:, 11:14] = A, D
    st[:, 14] = np.sqrt(np.einsum("ij,ij->i", D, D))
    st[:, 15] = np.arange(n)

    hints = getattr(shape, "_hints", None)
    with np.errstate(invalid="ignore", divide="ignore"):
        if hints is not None:
            # Segments the ends alone cannot rule out are measured next
            # wherever the shape expects their closest point, and cut there.
            keep = np.flatnonzero(_bound(st)[0] < min(best - tol, cutoff))
            if not keep.size:
                return best, who
            st = st[keep]
            H = hints(st[:, 8:11], st[:, 11:14])
            row, col = np.nonzero(np.isfinite(H) & (H > 0.0) & (H < 1.0))
            if row.size:
                th = H[row, col]
                fh, hh, gh = probe(st[row, 8:11] + th[:, None] * st[row, 11:14])
                sh = np.einsum("ij,ij->i", gh, st[row, 11:14])
                i = int(np.argmin(fh))
                if fh[i] < best:
                    best, who = float(fh[i]), int(st[row[i], 15])
                st = _cut(st, row, th, fh, hh, sh)

        for _ in range(_MAX_ROUNDS):
            lb, ta, tb = _bound(st)
            keep = np.flatnonzero(lb < min(best - tol, cutoff))
            if not keep.size:
                break
            if keep.size > _MAX_INTERVALS // 2:
                keep = keep[np.argsort(lb[keep])[:_MAX_INTERVALS // 2]]
            st = st[keep]
            # Two new points per stretch, so each becomes up to three.
            k = len(st)
            row = np.concatenate([np.arange(k), np.arange(k)])
            t = np.concatenate([ta[keep], tb[keep]])
            fm, hm, gm = probe(st[row, 8:11] + t[:, None] * st[row, 11:14])
            sm = np.einsum("ij,ij->i", gm, st[row, 11:14])
            i = int(np.argmin(fm))
            if fm[i] < best:
                best, who = float(fm[i]), int(st[row[i], 15])
            st = _cut(st, row, t, fm, hm, sm)
    return best, who


def _bound(st: np.ndarray):
    """Lower bound on the distance within each stretch, and two places to
    split it.

    See `_min_along` for the stretch layout and the two bounds.
    """
    t0, t1, f0, f1, h0, h1, s0, s1 = st[:, :8].T
    dt = t1 - t0
    lip = 0.5 * (f0 + f1 - dt * st[:, 14])
    # The convex bound is at least the larger of the tangent lines at the two
    # ends, which is lowest at an end or where they cross.
    # (np.minimum/np.maximum rather than np.clip, which costs several times
    # more per call at these sizes, and this runs every round.)
    x = np.minimum(np.maximum((h1 - h0 - s1 * dt) / (s0 - s1), 0.0), dt)
    x = np.where(np.isfinite(x), x, 0.0)
    hb = np.minimum(np.minimum(np.maximum(h0, h1 - s1 * dt),
                               np.maximum(h0 + s0 * dt, h1)),
                    np.maximum(h0 + s0 * x, h1 + s1 * (x - dt)))
    # Where to measure next, when the bound dips inside the stretch: two
    # guesses at the minimum, each exact for one kind of function.
    #   - where the tangents cross: exact where the distance is piecewise
    #     linear, as it is inside a solid, kinked where two faces are equally
    #     deep. For a smooth minimum it is useless — for a quadratic the
    #     tangents always cross at the midpoint — so on its own it only
    #     halves the stretch, a bit per round.
    #   - where the slope, interpolated linearly between the ends, is zero:
    #     exact for a quadratic, blind to a kink.
    # Both are kept off the ends so every split makes progress. When the
    # bound is lowest at an end, the stretch is alive only because the true
    # distance there exceeds the bound — the shadow of a bore — and halving
    # is the better guess.
    inner = (hb > lip) & (s0 < 0) & (s1 > 0)
    lo, hi = dt / 16, dt - dt / 16
    xs = -s0 * dt / (s1 - s0)
    ta = t0 + np.where(inner, np.minimum(np.maximum(x, lo), hi), 0.5 * dt)
    tb = t0 + np.where(inner, np.minimum(np.maximum(xs, lo), hi), 0.5 * dt)
    return np.maximum(lip, hb), ta, tb


def _cut(st: np.ndarray, row, t, f, h, s) -> np.ndarray:
    """Cut stretches at extra measured points: point k lies in stretch row[k].

    Every stretch's two ends and its extra points are sorted along it, and
    each consecutive pair becomes a stretch of its own.
    """
    n = len(st)
    idx = np.arange(n)
    R = np.concatenate([idx, idx, row])
    T = np.concatenate([st[:, 0], st[:, 1], t])
    F = np.concatenate([st[:, 2], st[:, 3], f])
    Hc = np.concatenate([st[:, 4], st[:, 5], h])
    S = np.concatenate([st[:, 6], st[:, 7], s])
    order = np.lexsort((T, R))
    R, T, F, Hc, S = R[order], T[order], F[order], Hc[order], S[order]
    a = np.flatnonzero(R[:-1] == R[1:])
    out = st[R[a]]
    out[:, 0], out[:, 2], out[:, 4], out[:, 6] = T[a], F[a], Hc[a], S[a]
    b = a + 1
    out[:, 1], out[:, 3], out[:, 5], out[:, 7] = T[b], F[b], Hc[b], S[b]
    return out


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
        # Long beads, by coarse cell. With two grids, a bead longer than
        # _LONG_BEAD_CELLS fine cells is filed in the coarse grid only, and
        # fine queries look here as well. See `add`.
        self._coarse_long: dict[tuple, list[int]] | None = (
            {} if self._coarse is not None else None)
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
        """Grid cells that any part of this bead reaches into.

        The bead's axis is sampled at most half a cell apart, and each sample
        padded by the bead radius plus half the spacing. Every point of the
        bead lies within `bead_radius` of some axis point, which lies within
        half a spacing of a sample, so no part of the bead falls outside the
        cells returned. That is linear in the bead's length.

        Two ways of doing this were wrong before it. Filing a bead under its
        whole bounding box is quadratic in its length: a 200 mm diagonal at
        2.75 mm cells is 5 000 cells, at needle-sized cells 150 000. Filing it
        under its two ends and midpoint only, which 0.4.0 did, matched a
        narrow phase that looked at nothing else, and inherited its blind
        spot: a needle in the middle of a long bead never saw it at all.
        """
        d = b - a
        length = math.sqrt(float(d @ d))
        n = max(1, math.ceil(2.0 * length / cell))
        pad = self.bead_radius + 0.5 * length / n
        pts = a + np.linspace(0.0, 1.0, n + 1)[:, None] * d
        lo = np.floor((pts - pad) / cell).astype(np.int64)
        hi = np.floor((pts + pad) / cell).astype(np.int64)
        if n <= 64:
            # Up to a few dozen samples, plain loops beat the array setup.
            out = set()
            for (l0, l1, l2), (h0, h1, h2) in zip(lo.tolist(), hi.tolist()):
                for i in range(l0, h0 + 1):
                    for j in range(l1, h1 + 1):
                        for k in range(l2, h2 + 1):
                            out.add((i, j, k))
            return out
        # Long bead: the same enumeration, vectorised. Each sample's box is
        # at most `span` cells a side, so offset every lower corner by every
        # step up to that, drop what overshoots the upper corner, and dedupe.
        span = int((hi - lo).max()) + 1
        steps = np.stack(np.meshgrid(*[np.arange(span)] * 3, indexing="ij"),
                         axis=-1).reshape(-1, 3)
        cells = (lo[:, None, :] + steps[None]).reshape(-1, 3)
        cells = cells[(cells <= np.repeat(hi, len(steps), axis=0)).all(axis=1)]
        base = cells.min(axis=0)
        rel = cells - base
        dims = rel.max(axis=0) + 1
        key = np.unique((rel[:, 0] * dims[1] + rel[:, 1]) * dims[2] + rel[:, 2])
        k2 = key % dims[2]
        k1 = (key // dims[2]) % dims[1]
        k0 = key // (dims[1] * dims[2])
        return set(zip((k0 + base[0]).tolist(), (k1 + base[1]).tolist(),
                       (k2 + base[2]).tolist()))

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
        # A bead is filed along its whole length, so its cost in cells is its
        # length over the cell size. At needle-sized cells that is what a
        # large part cannot afford: a 240 mm tray's infill lines put 1.5
        # million cells in the fine grid over 18 000 beads, and filing them
        # was half the run. So when there are two grids, a long bead goes in
        # the coarse one only — seven times fewer cells — and fine queries
        # also look up long beads by coarse cell. They get back more beads
        # than they would from fine cells, but the narrow phase discards a
        # distant bead on its first bound, which is cheap; filing it finely
        # is not.
        long = (self._coarse is not None
                and float(np.linalg.norm(B - A)) > _LONG_BEAD_CELLS * self.cell)
        if not long:
            for c in self._cells_for(A, B, self.cell):
                self._grid.setdefault(c, []).append(idx)
        if self._coarse is not None:
            cells = self._cells_for(A, B, self.coarse_cell)
            for c in cells:
                self._coarse.setdefault(c, []).append(idx)
            if long:
                for c in cells:
                    self._coarse_long.setdefault(c, []).append(idx)

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
        self.last_probe = self._collect(grid, lo_i, hi_i, n_cells, out)
        if grid is self._grid and self._coarse_long:
            # Long beads live in the coarse grid only; see `add`.
            lo_c, hi_c, n_c = self._box(lo, hi, self.coarse_cell)
            self.last_probe += self._collect(self._coarse_long, lo_c, hi_c,
                                             n_c, out)
        return np.fromiter(out, dtype=int, count=len(out))

    @staticmethod
    def _collect(grid, lo_i, hi_i, n_cells, out: set) -> int:
        """Add to `out` every index in `grid` cells within [lo_i, hi_i].

        Walks whichever is smaller, the box's cells or the occupied ones; both
        give the same set. Returns how many cells it inspected.
        """
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
        return min(n_cells, len(grid))

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
            search: broad-phase radius, metres. Every bead that comes within
                `search` of the shape is measured, along its whole length, so
                a result <= `search` is exact to `_NARROW_TOL`. A result above
                it is only whatever the grid happened to hand back — the
                nearest bead may not have been looked at — and is not a
                distance to rely on. `simulate` widens the search until it
                gets an exact answer, or reports `inf`.

        Returns:
            (clearance, segment_index). `inf` and -1 when nothing is in range.
            The clearance is a distance actually measured to a point on that
            bead, so it is never less than the true one.
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

        r = self.bead_radius
        d, j = _min_along(shape, self._A[cand], self._B[cand],
                          cutoff=search + r)
        return d - r, int(cand[j])

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
