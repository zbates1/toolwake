"""Tool profiles: a coaxial stack of sections, the way CAM describes a tool.

A dispensing needle is not a single cylinder. A luer-lock tip is a thin cannula
that steps up to a hub several tens of times its diameter, and the hub is what
actually fouls a tall part — exactly the holder-vs-flute distinction a CAM
simulation draws.

Sections are listed from the TIP upward, each with a length and a radius at
each end, so a straight cylinder and a taper are the same primitive.

For collision each section becomes a capsule at its LARGEST radius. That is
deliberately conservative: on a bioprinter a false stop costs a reprint and a
missed collision costs the part and possibly the needle. Tapers can be sliced
into sub-capsules when the loose fit matters more than the caution.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import Cylinder

__all__ = ["Section", "ToolProfile", "luer_taper_tip", "blunt_cannula"]


@dataclass(frozen=True)
class Section:
    """One coaxial piece of the tool.

    Args:
        length: axial length, metres.
        r0: radius at the lower (tip-ward) end, metres.
        r1: radius at the upper end, metres. Defaults to `r0` (a cylinder).
        name: label, used in reports so a collision says which part hit.
        r_inner: bore radius, metres. Non-zero makes the section a TUBE, which
            is what a cannula is. Material under an open bore is in the hole,
            not against the wall — model the needle solid and the layer it is
            printing onto reads as a collision on every row.
    """

    length: float
    r0: float
    r1: float | None = None
    name: str = ""
    r_inner: float = 0.0

    @property
    def r_top(self) -> float:
        return self.r0 if self.r1 is None else self.r1

    @property
    def r_max(self) -> float:
        return max(self.r0, self.r_top)

    def __post_init__(self):
        if self.length < 0:
            raise ValueError(f"section {self.name!r}: length must be >= 0")
        if self.r0 < 0 or self.r_top < 0:
            raise ValueError(f"section {self.name!r}: radii must be >= 0")
        if self.r_inner < 0:
            raise ValueError(f"section {self.name!r}: r_inner must be >= 0")
        if self.r_inner >= min(self.r0, self.r_top) and self.r_inner > 0:
            raise ValueError(
                f"section {self.name!r}: r_inner ({self.r_inner}) must be "
                f"smaller than the wall radius")


class ToolProfile:
    """An ordered stack of sections, tip first.

    The profile is expressed in the TOOL frame with +Z running back up the tool
    away from the tip, and placed into the world by `capsules()`.
    """

    def __init__(self, sections):
        self.sections = list(sections)
        if not self.sections:
            raise ValueError("a tool profile needs at least one section")

    # ------------------------------------------------------------- geometry
    @property
    def length(self) -> float:
        return float(sum(s.length for s in self.sections))

    @property
    def r_max(self) -> float:
        return float(max(s.r_max for s in self.sections))

    @property
    def tip_radius(self) -> float:
        """Radius right at the tip — what lays the bead."""
        return float(self.sections[0].r0)

    def stations(self) -> np.ndarray:
        """(M, 2) of (z, r) along the silhouette, from tip upward.

        A step between sections appears as two samples at the same z, so the
        outline drawn from this has square shoulders rather than ramps.
        """
        z, out = 0.0, [(0.0, self.sections[0].r0)]
        for s in self.sections:
            if abs(s.r0 - out[-1][1]) > 1e-12:      # step change in radius
                out.append((z, s.r0))
            z += s.length
            out.append((z, s.r_top))
        return np.asarray(out, dtype=float)

    def capsules(self, tip, axis, slices: int = 1):
        """World-space solids for collision, one (or `slices`) per section.

        Returns flat-ended `Cylinder`s, NOT capsules, despite the name (kept
        for compatibility). A capsule's hemispherical cap bulges a full radius
        past its endpoint, which for the tip-most section means the model
        occupies space below the nozzle that the needle does not — and the
        material directly below the nozzle is the print. See `Cylinder`.

        Args:
            tip: world position of the tool tip.
            axis: unit vector pointing from the tip back up the tool.
            slices: subdivisions per tapered section. 1 keeps each section at
                its largest radius (conservative); higher follows the taper
                more tightly at proportional cost.

        Returns:
            list of (name, Cylinder).
        """
        tip = np.asarray(tip, dtype=float)
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        out, z = [], 0.0
        for s in self.sections:
            if s.length <= 0:
                continue
            n = max(1, int(slices)) if s.r1 is not None and s.r1 != s.r0 else 1
            for k in range(n):
                t0, t1 = k / n, (k + 1) / n
                za, zb = z + s.length * t0, z + s.length * t1
                ra = s.r0 + (s.r_top - s.r0) * t0
                rb = s.r0 + (s.r_top - s.r0) * t1
                out.append((s.name or "section",
                            Cylinder(tip + axis * za, tip + axis * zb,
                                     max(ra, rb), s.r_inner)))
            z += s.length
        return out

    def surface(self, tip, R, n_theta: int = 24):
        """(X, Y, Z) arrays for drawing the tool as a surface of revolution."""
        st = self.stations()
        th = np.linspace(0.0, 2 * np.pi, n_theta)
        zz = st[:, 0][:, None]
        rr = st[:, 1][:, None]
        local = np.stack([rr * np.cos(th), rr * np.sin(th),
                          np.broadcast_to(zz, (len(st), n_theta))], axis=-1)
        world = local @ R.T + np.asarray(tip, dtype=float)
        return world[..., 0], world[..., 1], world[..., 2]

    def __repr__(self):
        names = ", ".join(s.name or "?" for s in self.sections)
        return (f"ToolProfile({len(self.sections)} sections: {names}; "
                f"len={self.length*1e3:.1f} mm, r_max={self.r_max*1e3:.2f} mm)")


# ------------------------------------------------------------------ presets
#
# NOMINAL dimensions. Hub sizes in particular vary between manufacturers and
# these are close enough to reason about clearance but not a substitute for
# calipers on the tip you are actually running.

def blunt_cannula(gauge_od=0.19e-3, length=12.7e-3, inner_d=90e-6) -> ToolProfile:
    """Just the metal tube — no hub. The optimistic model, for comparison.

    `inner_d` is the bore, and it is load-bearing: it was accepted and silently
    discarded here, so the cannula was modelled as a solid rod. A solid rod has
    no hole for the bead to come out of, and counts the layer it is printing
    onto as an obstacle.
    """
    return ToolProfile([Section(length, gauge_od / 2.0, name="cannula",
                                r_inner=inner_d / 2.0)])


def luer_taper_tip(cannula_od=0.19e-3, cannula_len=12.7e-3,
                   hub_d=9.0e-3, hub_len=13.0e-3,
                   taper_len=3.0e-3, collar_d=11.0e-3,
                   collar_len=4.0e-3, inner_d=90e-6) -> ToolProfile:
    """A luer-lock dispensing tip: cannula, taper, hub, locking collar.

    Defaults approximate a 34G half-inch tip. The hub is roughly 47x the
    cannula diameter, which is the entire reason this class exists — a check
    that models only the cannula is looking at the thinnest 2% of the tool.

    Args:
        cannula_od: needle outer diameter, metres (34G ~ 0.19 mm).
        cannula_len: exposed needle length, metres (0.5 in = 12.7 mm).
        hub_d: widest diameter of the plastic hub, metres.
        hub_len: straight hub length, metres.
        taper_len: cone joining cannula to hub, metres.
        collar_d, collar_len: the luer locking collar above the hub.
        inner_d: bore diameter, metres. Only the cannula is modelled hollow —
            above it the lumen is a rounding error against a 9 mm hub, and
            nothing that far up the tool is ever near a fresh bead.
    """
    return ToolProfile([
        Section(cannula_len, cannula_od / 2.0, name="cannula",
                r_inner=inner_d / 2.0),
        Section(taper_len, cannula_od / 2.0, hub_d / 2.0, name="taper"),
        Section(hub_len, hub_d / 2.0, name="hub"),
        Section(collar_len, collar_d / 2.0, name="luer collar"),
    ])
