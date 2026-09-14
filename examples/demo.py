"""Two runs: one that clears, one that does not.

The second is the interesting one — it is the case a needle-tip-only check
passes and a housing-aware check catches.
"""
import numpy as np

from toolwake import Needle, Toolpath, animate, simulate
from toolwake.toolpath import PRINT, TRAVEL


def clean():
    """A helix with enough pitch that the tool never revisits earlier material."""
    path = Toolpath.helix(radius=0.018, pitch=0.0035, turns=7, per_turn=110)
    needle = Needle(inner_d=90e-6, outer_d=0.4e-3, length=12.7e-3)
    return path, needle


def housing_collision():
    """Print a tall wall, then travel back down beside it.

    25 mm away: a 0.4 mm needle is comfortably clear, a 130 mm housing is not.
    """
    n = 150
    wall = np.column_stack([np.zeros(n), np.zeros(n), np.linspace(0, 0.06, n)])
    m = 60
    aside = np.column_stack([np.full(m, 0.025), np.zeros(m),
                             np.linspace(0.06, 0.0, m)])
    path = Toolpath(np.vstack([wall, aside]),
                    kinds=[PRINT] * n + [TRAVEL] * m)
    needle = Needle(inner_d=90e-6, outer_d=0.4e-3, length=12.7e-3,
                    housing=(0.130, 0.075, 0.072))
    return path, needle


if __name__ == "__main__":
    for name, build in (("helix_wake", clean),
                        ("housing_collision", housing_collision)):
        path, needle = build()
        res = simulate(path, needle, lag=10)
        print(f"{name}: {res!r}")
        print(f"  {res.report()}")
        animate(res, f"examples/{name}.mp4", fps=30, max_frames=300)


def nonplanar_two_walls():
    """Two concentric walls; the second is printed with the tool tilted outward.

    This is the realistic non-planar failure. A single vase wall rarely traps
    its own tool — the body leans over empty space. Put a second feature beside
    it and the leaning body sweeps straight through the one printed first.

    Outer wall at r = 20 mm goes down first. The inner wall at r = 10 mm is
    then printed with a 30 deg outward lean, which throws the hub about
    20*sin(30) = 10 mm outward — exactly onto the outer wall.
    """
    import numpy as np
    from toolwake import Needle, Toolpath
    from toolwake.geometry import rotvec_from_axis
    from toolwake.toolpath import PRINT, TRAVEL

    def ring(radius, z0, z1, turns, per_turn, lean):
        n = int(turns * per_turn)
        t = np.linspace(0.0, turns * 2 * np.pi, n)
        xyz = np.column_stack([radius * np.cos(t), radius * np.sin(t),
                               np.linspace(z0, z1, n)])
        axes = np.column_stack([np.sin(lean) * np.cos(t),
                                np.sin(lean) * np.sin(t),
                                np.full(n, np.cos(lean))])
        return xyz, np.array([rotvec_from_axis(u) for u in axes])

    outer_xyz, outer_rv = ring(0.020, 0.0, 0.030, 5, 80, 0.0)
    inner_xyz, inner_rv = ring(0.010, 0.0, 0.030, 5, 80, np.radians(30.0))

    # One travel hop between the two walls, which deposits nothing.
    hop_xyz = np.array([inner_xyz[0]])
    hop_rv = np.array([inner_rv[0]])

    path = Toolpath(np.vstack([outer_xyz, hop_xyz, inner_xyz]),
                    kinds=[PRINT] * len(outer_xyz) + [TRAVEL] + [PRINT] * len(inner_xyz),
                    rotvec=np.vstack([outer_rv, hop_rv, inner_rv]))
    return path, Needle.luer(inner_d=90e-6)
