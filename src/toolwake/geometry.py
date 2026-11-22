"""Exact distance primitives for the tool side of the check.

The tool is modelled with primitives rather than meshes on purpose. On the rig
this was written for, the end-effector's collision geometry in the URDF is
already a box and the needle is a thin cylinder, so a capsule and an oriented
box reproduce the real collision volume exactly — nothing is approximated away,
and there is no mesh library to depend on.

All functions are vectorised over a set of query points, because the caller
tests one tool pose against many deposited points at once.
"""
from __future__ import annotations

import numpy as np

__all__ = ["Capsule", "Cylinder", "Box", "rotation_from_rotvec", "rotvec_from_axis",
           "rpy_from_rotation"]


def rotation_from_rotvec(rotvec) -> np.ndarray:
    """Rotation matrix from a rotation vector (axis * angle), via Rodrigues.

    Matches the convention UR controllers use for pose orientation, so a pose
    row read straight off the robot can be passed in unchanged.
    """
    r = np.asarray(rotvec, dtype=float)
    theta = float(np.linalg.norm(r))
    if theta < 1e-12:
        return np.eye(3)
    k = r / theta
    K = np.array([[0.0, -k[2], k[1]],
                  [k[2], 0.0, -k[0]],
                  [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def rotvec_from_axis(u, reference=(0.0, 0.0, 1.0)) -> np.ndarray:
    """Rotation vector taking `reference` onto `u` by the shortest arc.

    The inverse of what `Needle.at` needs: give it the direction the tool
    should point (tip-to-body, so "up the tool") and it returns the rotation
    vector to store on the toolpath.

    The shortest arc leaves the roll about the tool axis unconstrained, which
    is correct here — a round needle laying a round bead has no meaningful
    roll, so pinning one would be inventing a constraint.
    """
    a = np.asarray(reference, dtype=float)
    b = np.asarray(u, dtype=float)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-15 or nb < 1e-15:
        raise ValueError("rotvec_from_axis needs two non-zero vectors")
    a, b = a / na, b / nb
    c = np.cross(a, b)
    s = float(np.linalg.norm(c))
    d = float(np.dot(a, b))
    if s < 1e-12:
        if d > 0:
            return np.zeros(3)
        # Antiparallel: any perpendicular axis gives a 180 deg flip.
        perp = np.array([1.0, 0.0, 0.0])
        if abs(a[0]) > 0.9:
            perp = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perp)
        return axis / np.linalg.norm(axis) * np.pi
    return c / s * float(np.arctan2(s, d))


def rpy_from_rotation(R) -> np.ndarray:
    """(roll, pitch, yaw) in radians from a rotation matrix, ZYX convention.

    Yaw about Z, then pitch about the new Y, then roll about the new X — the
    convention robot controllers and teach pendants report, so the numbers here
    match what an operator reads off the machine.

    At pitch = +/-90 deg roll and yaw describe the same rotation (gimbal lock);
    there the split is arbitrary and roll is pinned to 0 so the result stays
    deterministic rather than amplifying numerical noise.
    """
    R = np.asarray(R, dtype=float)
    sp = -R[2, 0]
    if abs(sp) > 1.0 - 1e-9:                     # gimbal lock
        pitch = np.pi / 2 * np.sign(sp)
        return np.array([0.0, pitch, float(np.arctan2(-R[0, 1], R[1, 1]))])
    return np.array([float(np.arctan2(R[2, 1], R[2, 2])),
                     float(np.arcsin(np.clip(sp, -1.0, 1.0))),
                     float(np.arctan2(R[1, 0], R[0, 0]))])


def _unit_rows(v: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Each row of `v` divided by its norm `n`; zero where the norm is zero.

    Zero is the right answer there, not an arbitrary one: these are gradients
    of convex functions, and at a point where the vector vanishes the function
    is at its minimum, where zero is a valid subgradient.
    """
    out = np.zeros_like(v)
    ok = n > 0
    out[ok] = v[ok] / n[ok, None]
    return out


def _seg_point_distance(pts: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Distance from each point to the SEGMENT ab (not the infinite line)."""
    pts = np.atleast_2d(pts)
    ab = b - a
    denom = float(ab @ ab)
    if denom < 1e-18:                      # degenerate segment -> a point
        return np.linalg.norm(pts - a, axis=1)
    t = np.clip((pts - a) @ ab / denom, 0.0, 1.0)
    closest = a + t[:, None] * ab
    return np.linalg.norm(pts - closest, axis=1)


class Capsule:
    """A cylinder with hemispherical caps: the needle, or a deposited bead.

    Defined by its axis segment and a radius. A capsule is used rather than a
    plain cylinder because the distance function stays exact and smooth at the
    ends, which matters when a bead is one short segment.
    """

    __slots__ = ("a", "b", "radius")

    def __init__(self, a, b, radius: float):
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.radius = float(radius)
        if self.a.shape != (3,) or self.b.shape != (3,):
            raise ValueError("capsule endpoints must be 3-vectors")
        if self.radius < 0:
            raise ValueError("capsule radius must be >= 0")

    def distance(self, pts) -> np.ndarray:
        """Signed distance from points to the capsule surface (negative = inside)."""
        return _seg_point_distance(pts, self.a, self.b) - self.radius

    def _probe(self, pts):
        """Distance, a convex lower bound on it, and that bound's gradient.

        See `Cylinder._probe` for what these are for. A capsule is convex, so
        its own distance is the bound.
        """
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        ab = self.b - self.a
        denom = float(ab @ ab)
        if denom < 1e-18:
            diff = pts - self.a
        else:
            t = np.clip((pts - self.a) @ ab / denom, 0.0, 1.0)
            diff = pts - (self.a + t[:, None] * ab)
        n = np.linalg.norm(diff, axis=1)
        f = n - self.radius
        return f, f, _unit_rows(diff, n)

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.b - self.a))

    def bounds(self, pad: float = 0.0):
        lo = np.minimum(self.a, self.b) - self.radius - pad
        hi = np.maximum(self.a, self.b) + self.radius + pad
        return lo, hi

    def __repr__(self):
        return f"Capsule(len={self.length:.4g}, r={self.radius:.4g})"


class Cylinder:
    """A cylinder with FLAT ends, optionally hollow: a tool section.

    The difference from `Capsule` is the caps, and it is not cosmetic. A
    capsule's hemispherical end bulges a full radius beyond its axis endpoint,
    so a cannula modelled as one reaches `radius` BELOW the nozzle tip — into
    space the needle does not occupy.

    That produced a specific, confident wrong answer. On a 0.06 mm-layer print
    with a 0.095 mm-radius cannula, the phantom hemisphere swallowed the two
    layers directly underneath the nozzle on essentially every row: 2 193 of
    2 674 rows reported as collisions, worst -0.14 mm, all blamed on the
    cannula. With flat caps the same print's tightest clearance is +0.015 mm.
    Material under the nozzle is what a print IS; a tool model must not treat
    it as an obstacle.

    `inner_radius` makes it a TUBE, which is what a cannula is. The bore is a
    hole, so material sitting under it — the bead the nozzle is laying, and the
    layer it is laying onto — is inside the hole, not touching the wall. Model
    the needle as solid and every one of those reads as a collision.

    Beads stay capsules — a short segment of extruded material really does have
    rounded ends.
    """

    __slots__ = ("a", "b", "radius", "inner_radius", "_L", "_u")

    def __init__(self, a, b, radius: float, inner_radius: float = 0.0):
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.radius = float(radius)
        self.inner_radius = float(inner_radius)
        if self.a.shape != (3,) or self.b.shape != (3,):
            raise ValueError("cylinder endpoints must be 3-vectors")
        # Axis length and direction, fixed at construction: the narrow phase
        # asks for them on every probe.
        self._L = float(np.linalg.norm(self.b - self.a))
        self._u = (self.b - self.a) / self._L if self._L >= 1e-18 else None
        if self.radius < 0:
            raise ValueError("cylinder radius must be >= 0")
        if self.inner_radius < 0:
            raise ValueError("cylinder inner_radius must be >= 0")
        if self.inner_radius >= self.radius > 0:
            raise ValueError("cylinder inner_radius must be < radius")

    def distance(self, pts) -> np.ndarray:
        """Signed distance to the surface (negative inside), exact.

        In (radial, axial) coordinates the solid is the rectangle
        [inner_radius, radius] x [0, L], so this is the standard 2-D rectangle
        SDF: outside, the hypotenuse of however far the point lies past the
        wall and past the end; inside, the distance to whichever surface is
        nearest, which is the LARGER (least negative) of the two.

        A SOLID cylinder has no inner surface, so it must not get one at the
        axis. Taking the rectangle literally as [0, radius] did exactly that:
        a point on the axis, deep inside a hub, measured as clearance 0 —
        touching, not penetrating — and depth near the axis read as the
        distance to the axis rather than to the wall.
        """
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        geo = self._local(pts)
        if geo is None:                     # degenerate -> a disc; treat as sphere
            return np.linalg.norm(pts - self.a, axis=1) - self.radius
        _, _, radial, dz = geo
        dr = radial - self.radius
        if self.inner_radius > 0:
            dr = np.maximum(self.inner_radius - radial, dr)
        outside = np.sqrt(np.maximum(dr, 0.0) ** 2 + np.maximum(dz, 0.0) ** 2)
        return np.where((dr <= 0) & (dz <= 0), np.maximum(dr, dz), outside)

    def _hints(self, A, D):
        """Where along each segment A + t D its closest point is likely, (n, 7).

        Where it passes nearest the axis, where it crosses the bore and outer
        radii, and where it crosses the two end planes. NaN where a crossing
        does not exist. These only seed the narrow phase's search; nothing is
        concluded from them, so a poor guess costs time, not correctness.

        The bore crossings are the ones that matter most. A bead passing under
        the needle is closest where it enters the wall's shadow, but the solid
        cylinder that bounds the tube is lowest on the axis, in the bore's
        shadow — so without these the search halves its way outward from the
        axis to the wall, seven rounds on a typical infill line.
        """
        n = len(A)
        if self._u is None:
            return np.full((n, 7), np.nan)
        u = self._u
        rel = A - self.a
        z0, dz = rel @ u, D @ u
        w0 = rel - z0[:, None] * u
        wd = D - dz[:, None] * u
        al = np.einsum("ij,ij->i", wd, wd)
        be = np.einsum("ij,ij->i", w0, wd)
        ga = np.einsum("ij,ij->i", w0, w0)
        out = np.full((n, 7), np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            out[:, 0] = -be / al
            for c, col in ((self.inner_radius, 1), (self.radius, 3)):
                root = np.sqrt(be * be - al * (ga - c * c))
                out[:, col] = (-be - root) / al
                out[:, col + 1] = (-be + root) / al
            out[:, 5] = -z0 / dz
            out[:, 6] = (self._L - z0) / dz
        return out

    def _local(self, pts):
        """(axial, radial vector, radial distance, axial excess), or None for
        a degenerate cylinder."""
        u = self._u
        if u is None:
            return None
        rel = pts - self.a
        t = rel @ u                                     # axial, 0..L is inside
        rvec = rel - t[:, None] * u
        radial = np.sqrt(np.einsum("ij,ij->i", rvec, rvec))
        dz = np.maximum(-t, t - self._L)
        return t, rvec, radial, dz

    def _probe(self, pts):
        """Distance, a convex lower bound on it, and that bound's gradient.

        `Deposit` uses these to find the closest point along a whole bead
        rather than at a few points on it. Distance alone only bounds what
        lies between two measured points by how far apart they are, and a
        bead lying flat under the nozzle face is at the same distance all
        along — so proving nothing on it is closer would mean measuring it
        every nanometre. A convex function bounds itself much more tightly:
        it never dips below the tangent lines at the two ends, so one gradient
        per point settles a flat stretch at once.

        A tube is not convex — the bore is a hole — so its bound is the SOLID
        cylinder it is drilled from. That never exceeds the tube's distance
        and equals it everywhere outside the bore's shadow, which is where the
        closest bead almost always is.
        """
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        geo = self._local(pts)
        if geo is None:
            diff = pts - self.a
            n = np.linalg.norm(diff, axis=1)
            f = n - self.radius
            return f, f, _unit_rows(diff, n)
        t, rvec, radial, dz = geo
        drs = radial - self.radius
        out_r, out_z = np.maximum(drs, 0.0), np.maximum(dz, 0.0)
        outside = np.sqrt(out_r * out_r + out_z * out_z)
        inside = (drs <= 0) & (dz <= 0)
        h = np.where(inside, np.maximum(drs, dz), outside)
        if self.inner_radius > 0:
            dr = np.maximum(self.inner_radius - radial, drs)
            out_b = np.maximum(dr, 0.0)
            f = np.where((dr <= 0) & (dz <= 0), np.maximum(dr, dz),
                         np.sqrt(out_b * out_b + out_z * out_z))
        else:
            f = h
        # Gradient of h: radial and axial weights, each a unit direction.
        # Axially, dz grows away from whichever end is nearer. On the axis
        # the radial direction is undefined and its weight is zero there —
        # h is at its radial minimum, where zero is a valid subgradient.
        safe = np.where(outside > 0, outside, 1.0)
        w_r = np.where(inside, drs > dz, out_r / safe)
        w_z = np.where(inside, drs <= dz, out_z / safe)
        w_z = np.where(t < 0.5 * self._L, -w_z, w_z)
        w_r = np.divide(w_r, radial, out=np.zeros_like(radial), where=radial > 0)
        return f, h, w_r[:, None] * rvec + w_z[:, None] * self._u

    @property
    def length(self) -> float:
        return float(np.linalg.norm(self.b - self.a))

    def bounds(self, pad: float = 0.0):
        lo = np.minimum(self.a, self.b) - self.radius - pad
        hi = np.maximum(self.a, self.b) + self.radius + pad
        return lo, hi

    def __repr__(self):
        bore = (f", inner_radius={self.inner_radius:.4g}"
                if self.inner_radius else "")
        return f"Cylinder(a={self.a}, b={self.b}, radius={self.radius:.4g}{bore})"


class Box:
    """An oriented box — the end-effector housing.

    `R` maps BOX-local coordinates into world. Distance is exact for points
    outside; for points inside it returns a negative penetration depth (the
    distance to the nearest face), which is the usual convention and is enough
    to rank how bad a collision is.
    """

    __slots__ = ("center", "half", "R")

    def __init__(self, center, half_extents, R=None):
        self.center = np.asarray(center, dtype=float)
        self.half = np.asarray(half_extents, dtype=float)
        self.R = np.eye(3) if R is None else np.asarray(R, dtype=float)
        if self.center.shape != (3,) or self.half.shape != (3,):
            raise ValueError("box center and half_extents must be 3-vectors")
        if np.any(self.half < 0):
            raise ValueError("box half_extents must be >= 0")

    def distance(self, pts) -> np.ndarray:
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        local = (pts - self.center) @ self.R          # world -> box frame
        excess = np.abs(local) - self.half
        outside = np.linalg.norm(np.maximum(excess, 0.0), axis=1)
        # Inside on every axis -> the closest face is the least-negative excess.
        inside = np.minimum(excess.max(axis=1), 0.0)
        return outside + inside

    def _probe(self, pts):
        """Distance, a convex lower bound on it, and that bound's gradient.

        See `Cylinder._probe`. A box is convex, so its distance is the bound.
        """
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        local = (pts - self.center) @ self.R
        excess = np.abs(local) - self.half
        pos = np.maximum(excess, 0.0)
        outside = np.linalg.norm(pos, axis=1)
        f = outside + np.minimum(excess.max(axis=1), 0.0)
        sgn = np.sign(local)
        # Outside: along the excess. Inside: out through the nearest face.
        nearest = (np.arange(3) == excess.argmax(axis=1)[:, None]).astype(float)
        g = np.where((outside > 0)[:, None], _unit_rows(sgn * pos, outside),
                     sgn * nearest)
        return f, f, g @ self.R.T

    def corners(self) -> np.ndarray:
        """The 8 world-space corners, for drawing a wireframe."""
        signs = np.array([[sx, sy, sz]
                          for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)],
                         dtype=float)
        return self.center + (signs * self.half) @ self.R.T

    def bounds(self, pad: float = 0.0):
        c = self.corners()
        return c.min(axis=0) - pad, c.max(axis=0) + pad

    def __repr__(self):
        return f"Box(half={np.round(self.half, 4).tolist()})"
