"""Toolpaths: where the tool goes, and which moves actually deposit.

The `kinds` array is the part that is easy to get wrong and expensive to get
wrong. If travel moves are treated as deposit, the wake fills with material the
printer never laid and every later move reports a phantom collision. So a
toolpath always carries a per-row kind, and only `print` rows leave a wake.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .geometry import rotvec_from_axis

__all__ = ["Toolpath", "PRINT", "TRAVEL"]

PRINT = "print"
TRAVEL = "travel"

_NUM = r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)"
_AXIS = {ax: re.compile(ax + _NUM) for ax in "XYZEF"}


class Toolpath:
    """Ordered tool positions plus, optionally, tool orientation.

    Attributes:
        xyz:   (N, 3) positions, metres.
        kinds: (N,) of "print" / "travel". Row i describes the move INTO xyz[i].
        rotvec: (N, 3) rotation vectors, or None for a planar path where the
            tool is assumed to point along -Z.
        feed:  (N,) feedrate in m/s, or None.
    """

    def __init__(self, xyz, kinds=None, rotvec=None, feed=None):
        self.xyz = np.asarray(xyz, dtype=float)
        if self.xyz.ndim != 2 or self.xyz.shape[1] != 3:
            raise ValueError(f"xyz must be (N, 3); got {self.xyz.shape}")
        n = len(self.xyz)
        self.kinds = (np.full(n, PRINT, dtype=object) if kinds is None
                      else np.asarray(kinds, dtype=object))
        if len(self.kinds) != n:
            raise ValueError("kinds must match xyz length")
        self.rotvec = None if rotvec is None else np.asarray(rotvec, dtype=float)
        if self.rotvec is not None and self.rotvec.shape != (n, 3):
            raise ValueError(f"rotvec must be (N, 3); got {self.rotvec.shape}")
        self.feed = None if feed is None else np.asarray(feed, dtype=float)

    def __len__(self):
        return len(self.xyz)

    @property
    def is_planar(self) -> bool:
        return self.rotvec is None

    def bounds(self):
        return self.xyz.min(axis=0), self.xyz.max(axis=0)

    def segment_lengths(self) -> np.ndarray:
        d = np.zeros(len(self))
        if len(self) > 1:
            d[1:] = np.linalg.norm(np.diff(self.xyz, axis=0), axis=1)
        return d

    # ------------------------------------------------------------- loaders
    @classmethod
    def from_arrays(cls, xyz, kinds=None, rotvec=None, feed=None) -> "Toolpath":
        return cls(xyz, kinds=kinds, rotvec=rotvec, feed=feed)

    @classmethod
    def from_gcode(cls, path, units: float = 1e-3) -> "Toolpath":
        """Parse G0/G1 moves from a G-code file.

        Deliberately minimal — enough to know position, whether the move
        extrudes, and its feedrate. `units` converts file units to metres
        (default: the file is in millimetres).

        Handles M82/M83 (absolute/relative E) because the distinction decides
        whether a row extrudes, and slicers disagree: Cura emits M82,
        PrusaSlicer M83.
        """
        text = Path(path).read_text(errors="replace").splitlines()
        pos = np.zeros(3)
        e_abs, e_rel = 0.0, False
        feed = 0.0
        xyz, kinds, feeds = [], [], []

        for raw in text:
            line = raw.split(";", 1)[0].strip()
            if not line:
                continue
            head = line.split()[0].upper()
            if head == "M82":
                e_rel = False
                continue
            if head == "M83":
                e_rel = True
                continue
            if head == "G92":
                m = _AXIS["E"].search(line)
                if m:
                    e_abs = float(m.group(1))
                continue
            if head not in ("G0", "G1", "G00", "G01"):
                continue

            new = pos.copy()
            for i, ax in enumerate("XYZ"):
                m = _AXIS[ax].search(line)
                if m:
                    new[i] = float(m.group(1)) * units
            m = _AXIS["F"].search(line)
            if m:
                feed = float(m.group(1)) * units / 60.0     # mm/min -> m/s

            de = 0.0
            m = _AXIS["E"].search(line)
            if m:
                v = float(m.group(1))
                de, e_abs = (v, e_abs + v) if e_rel else (v - e_abs, v)

            moved = float(np.linalg.norm(new - pos)) > 0
            if moved:
                xyz.append(new.copy())
                kinds.append(PRINT if de > 0 else TRAVEL)
                feeds.append(feed)
            pos = new

        if not xyz:
            raise ValueError(f"no G0/G1 moves found in {path}")
        return cls(np.asarray(xyz), kinds=kinds, feed=np.asarray(feeds))

    @classmethod
    def from_poses(cls, path, units: float = 1.0, kinds=None) -> "Toolpath":
        """Load an (N, 6) pose file: x y z rx ry rz, one row per waypoint.

        Whitespace, comma or tab separated; `#` starts a comment. This is the
        shape non-planar pipelines emit and the shape a UR reports, so rows can
        move between the robot and here without conversion.

        Args:
            units: multiplier on the POSITION columns only. Rotation vectors
                are radians either way, so scaling them would be wrong.
            kinds: per-row print/travel. Every row is a print move if omitted,
                which is the usual case for a generated non-planar path.
        """
        rows = []
        for raw in Path(path).read_text(errors="replace").splitlines():
            line = raw.split("#", 1)[0].replace(",", " ").strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                rows.append([float(v) for v in parts[:6]])
            except ValueError:
                continue                    # header or stray text
        if not rows:
            raise ValueError(f"no 6-column pose rows found in {path}")
        arr = np.asarray(rows, dtype=float)
        return cls(arr[:, :3] * units, kinds=kinds, rotvec=arr[:, 3:6])

    def with_tool_axis(self, axes) -> "Toolpath":
        """Return a copy whose tool orientation follows `axes`.

        `axes` is (N, 3) of directions pointing FROM the tip back up the tool —
        the surface normal, for a path that should stay perpendicular to what
        it is printing on. This is how a planar path becomes a non-planar one.
        """
        axes = np.asarray(axes, dtype=float)
        if axes.shape != (len(self), 3):
            raise ValueError(f"axes must be ({len(self)}, 3); got {axes.shape}")
        rv = np.array([rotvec_from_axis(u) for u in axes])
        return Toolpath(self.xyz, kinds=self.kinds, rotvec=rv, feed=self.feed)

    @classmethod
    def conical_helix(cls, r0=0.008, r1=0.018, height=0.030, turns=8,
                      per_turn=90, tilt=None, z0=0.0) -> "Toolpath":
        """A flaring vase wall — the standard non-planar test case.

        The radius ramps from `r0` to `r1` over `height`, so the wall leans out
        and the tool has to lean with it. `tilt` is the angle from vertical, in
        radians; the default is the wall's own slope, which keeps the needle
        perpendicular to the surface it is printing on.

        Set `tilt=0` for the same geometry printed with a vertical tool — the
        two runs side by side are what makes non-planar support visible, since
        a leaning tool reaches over material a vertical one never approaches.
        """
        n = int(turns * per_turn)
        t = np.linspace(0.0, turns * 2 * np.pi, n)
        frac = t / (turns * 2 * np.pi)
        r = r0 + (r1 - r0) * frac
        z = z0 + height * frac
        xyz = np.column_stack([r * np.cos(t), r * np.sin(t), z])

        slope = (r1 - r0) / height if height > 0 else 0.0
        lean = float(np.arctan(slope)) if tilt is None else float(tilt)
        # Tool axis: vertical, leaned outward in each row's own radial plane.
        axes = np.column_stack([np.sin(lean) * np.cos(t),
                                np.sin(lean) * np.sin(t),
                                np.full(n, np.cos(lean))])
        return cls(xyz, kinds=[PRINT] * n).with_tool_axis(axes)

    @classmethod
    def helix(cls, radius=0.02, pitch=0.004, turns=6.0, per_turn=120,
              z0=0.0) -> "Toolpath":
        """A helical non-planar path — a compact stand-in for a vessel wall.

        Useful as a demo and as a test fixture: it is the simplest path whose
        later passes come back around directly above earlier ones, which is
        exactly the geometry a wake check has to reason about.
        """
        n = int(turns * per_turn)
        t = np.linspace(0.0, turns * 2 * np.pi, n)
        xyz = np.column_stack([radius * np.cos(t),
                               radius * np.sin(t),
                               z0 + pitch * t / (2 * np.pi)])
        return cls(xyz, kinds=[PRINT] * n)

    def __repr__(self):
        n_print = int(np.sum(self.kinds == PRINT))
        kind = "planar" if self.is_planar else "6-DOF"
        return (f"Toolpath({len(self)} rows, {n_print} print, {kind}, "
                f"{self.segment_lengths().sum() * 1e3:.0f} mm)")
