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

    __slots__ = ("a", "b", "radius", "inner_radius")

    def __init__(self, a, b, radius: float, inner_radius: float = 0.0):
        self.a = np.asarray(a, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.radius = float(radius)
        self.inner_radius = float(inner_radius)
        if self.a.shape != (3,) or self.b.shape != (3,):
            raise ValueError("cylinder endpoints must be 3-vectors")
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
        """
        pts = np.atleast_2d(np.asarray(pts, dtype=float))
        ab = self.b - self.a
        L = float(np.linalg.norm(ab))
        if L < 1e-18:                       # degenerate -> a disc; treat as sphere
            return np.linalg.norm(pts - self.a, axis=1) - self.radius
        u = ab / L
        rel = pts - self.a
        t = rel @ u                                     # axial, 0..L is inside
        radial = np.linalg.norm(rel - t[:, None] * u, axis=1)
        dr = np.maximum(self.inner_radius - radial, radial - self.radius)
        dz = np.maximum(-t, t - L)
        outside = np.sqrt(np.maximum(dr, 0.0) ** 2 + np.maximum(dz, 0.0) ** 2)
        return np.where((dr <= 0) & (dz <= 0), np.maximum(dr, dz), outside)

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
