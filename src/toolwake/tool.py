"""The tool: a profile of coaxial sections, optionally under a housing box.

Why the shape matters. The needle tip crossing an earlier bead is the failure
people picture, and it is the rarest one. A luer tip's hub is roughly 47x the
cannula diameter, and the end-effector body is larger again — on a tall part
those are inside the structure while the tip is somewhere harmless. A check
that models only the tip passes moves that physically cannot be made.

Geometry is expressed in the TOOL frame with +Z running back up the tool away
from the tip, and placed into the world by `at()`.
"""
from __future__ import annotations

import numpy as np

from .geometry import Box, Capsule, rotation_from_rotvec
from .profile import Section, ToolProfile, blunt_cannula, luer_taper_tip

__all__ = ["Needle", "ToolPose"]

_DOWN = np.array([0.0, 0.0, -1.0])


class ToolPose:
    """A placed tool: named world-space shapes, ready to test."""

    __slots__ = ("parts", "housing", "tip", "R", "axis")

    def __init__(self, parts, housing, tip, R, axis):
        self.parts = parts            # list of (name, Capsule), tip first
        self.housing = housing        # Box or None
        self.tip = tip
        self.R = R
        self.axis = axis

    @property
    def needle(self) -> Capsule:
        """The tip-most section — kept so older call sites still work."""
        return self.parts[0][1]

    def named_shapes(self):
        """Every testable shape with its label, housing last."""
        out = list(self.parts)
        if self.housing is not None:
            out.append(("housing", self.housing))
        return out

    def shapes(self):
        return [s for _, s in self.named_shapes()]


class Needle:
    """A dispensing tool.

    Two ways to build one:

        Needle(inner_d=90e-6, outer_d=0.4e-3, length=12.7e-3)   # plain cannula
        Needle.luer(inner_d=90e-6)                              # full luer tip

    Args:
        inner_d: needle inner diameter, metres — this sets the bead laid down.
        outer_d: cannula outer diameter, metres. Defaults to 2x inner.
        length: exposed cannula length, metres.
        housing: (x, y, z) full extents of the end-effector body, metres, or
            None to model the tool alone.
        housing_offset: centre of the housing in the tool frame; defaults to
            sitting directly above the tool profile.
        profile: an explicit `ToolProfile`, overriding outer_d/length.
    """

    def __init__(self, inner_d=90e-6, outer_d=None, length=12.7e-3,
                 housing=None, housing_offset=None, profile=None):
        self.inner_d = float(inner_d)
        self.outer_d = float(outer_d) if outer_d else 2.0 * self.inner_d
        self.profile = profile or blunt_cannula(
            gauge_od=self.outer_d, length=float(length), inner_d=self.inner_d)
        self.housing = None if housing is None else np.asarray(housing, float)
        if self.housing is not None and housing_offset is None:
            housing_offset = (0.0, 0.0, self.profile.length + self.housing[2] / 2.0)
        self.housing_offset = (None if housing_offset is None
                               else np.asarray(housing_offset, float))

    # ------------------------------------------------------------ builders
    @classmethod
    def luer(cls, inner_d=90e-6, housing=None, **kw) -> "Needle":
        """A luer-lock dispensing tip: cannula, taper, hub, collar.

        Extra keywords go to `profile.luer_taper_tip` — pass measured hub
        dimensions there rather than trusting the nominal defaults.
        """
        prof = luer_taper_tip(**kw)
        return cls(inner_d=inner_d, outer_d=prof.tip_radius * 2.0,
                   profile=prof, housing=housing)

    # ------------------------------------------------------------ geometry
    @property
    def length(self) -> float:
        return self.profile.length

    @property
    def bead_radius(self) -> float:
        """Radius of the material left behind — set by the INNER diameter."""
        return self.inner_d / 2.0

    def at(self, tip, rotvec=None, slices: int = 1) -> ToolPose:
        """Place the tool with its tip at `tip`.

        With no `rotvec` the tool points straight down (the planar case). With
        one, it is read the way a UR controller reports pose orientation, so a
        row taken off the robot can be passed in unchanged.
        """
        tip = np.asarray(tip, dtype=float)
        R = np.eye(3) if rotvec is None else rotation_from_rotvec(rotvec)
        axis = R @ -_DOWN                    # tool +Z, back up away from the tip
        parts = self.profile.capsules(tip, axis, slices=slices)
        housing = None
        if self.housing is not None:
            housing = Box(tip + R @ self.housing_offset, self.housing / 2.0, R)
        return ToolPose(parts, housing, tip, R, axis)

    def __repr__(self):
        h = "none" if self.housing is None else np.round(self.housing * 1e3, 1).tolist()
        return (f"Needle(id={self.inner_d*1e6:.0f}um, "
                f"profile={self.profile!r}, housing_mm={h})")
