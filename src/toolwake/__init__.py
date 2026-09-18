"""toolwake — simulate a depositing tool, accumulate its wake, check against it.

    from toolwake import Toolpath, Needle, simulate, animate

    path = Toolpath.helix(radius=0.02, pitch=0.004, turns=6)
    res  = simulate(path, Needle(inner_d=90e-6, housing=(0.130, 0.075, 0.072)))
    print(res.report())
    animate(res, "wake.mp4")

The wake the video draws and the wake the clearance numbers are measured
against are the same object, so the picture cannot disagree with the report.
"""
from .deposit import Deposit
from .geometry import (Box, Capsule, rotation_from_rotvec,
                       rotvec_from_axis)
from .simulate import Result, simulate
from .viewer import to_html_str
from .profile import Section, ToolProfile, blunt_cannula, luer_taper_tip
from .tool import Needle, ToolPose
from .toolpath import PRINT, TRAVEL, Toolpath

__version__ = "0.1.1"

__all__ = [
    "Toolpath", "PRINT", "TRAVEL",
    "Needle", "ToolPose",
    "Section", "ToolProfile", "blunt_cannula", "luer_taper_tip",
    "Deposit",
    "Capsule", "Box", "rotation_from_rotvec", "rotvec_from_axis",
    "simulate", "Result",
    "animate", "to_html", "to_html_str", "ToolwakeSession",
    "__version__",
]


def animate(*args, **kwargs):
    """Render a Result to video. Imported lazily so matplotlib stays optional."""
    from .render import animate as _animate
    return _animate(*args, **kwargs)


def ToolwakeSession(*args, **kwargs):
    """A re-checkable simulation served as HTML, for embedding in another GUI."""
    from .serve import ToolwakeSession as _S
    return _S(*args, **kwargs)


def to_html(*args, **kwargs):
    """Write a standalone interactive viewer. No dependencies beyond numpy."""
    from .viewer import to_html as _to_html
    return _to_html(*args, **kwargs)
