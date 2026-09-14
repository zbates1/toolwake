"""Animate the tool travelling its path and leaving a wake behind it.

Rendering is deliberately a view onto the same `Result` the collision check
produced — not a second model of the process. If the video shows the needle
passing through the wake, the numbers say so too.

matplotlib only, so the package stays installable without a GPU stack.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = ["animate"]

_NEEDLE = "#2D2D8A"
_WAKE = "#3FB8D4"
_HOT = "#C41230"
_TIP = "#FF9100"
_GRID = "#D8DDE4"


def _equal_box(ax, lo, hi, pad=0.08):
    """Equal aspect on all three axes — otherwise a helix renders as an ellipse."""
    span = float(np.max(hi - lo))
    span = span if span > 0 else 1.0
    c = 0.5 * (lo + hi)
    h = 0.5 * span * (1.0 + pad)
    ax.set_xlim(c[0] - h, c[0] + h)
    ax.set_ylim(c[1] - h, c[1] + h)
    ax.set_zlim(c[2] - h, c[2] + h)
    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect((1, 1, 1))


def animate(result, out="wake.mp4", *, fps=30, max_frames=600, dpi=140,
            figsize=(9, 6), elev=22.0, azim0=-60.0, spin=25.0,
            show_housing=True, title=None, progress=None, n_theta=20):
    """Render `result` to a video file.

    Args:
        result: a `Result` from `toolwake.simulate`.
        out: output path. `.mp4` needs ffmpeg; `.gif` uses Pillow.
        fps: frames per second of the output.
        max_frames: cap on rendered frames — the path is subsampled to fit, so
            a 40k-row toolpath still renders in a sensible time.
        spin: degrees of azimuth rotation across the whole clip. 0 to hold still.
        show_housing: draw the end-effector box when the tool has one.
        n_theta: facets around the tool's surface of revolution. Lower is
            faster; 20 is smooth enough at video resolution.
        progress: optional callable(i, n).

    Returns:
        The path written.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FFMpegWriter, FuncAnimation, PillowWriter
    from mpl_toolkits.mplot3d.art3d import Line3DCollection

    path, needle = result.path, result.needle
    n = len(path)
    step = max(1, int(np.ceil(n / max_frames)))
    frames = list(range(0, n, step))
    if frames[-1] != n - 1:
        frames.append(n - 1)

    lo, hi = path.bounds()
    # Pad for the TOOL, not just the path: a 130 mm housing drawn around a
    # 60 mm part is otherwise sliced off by the axes and the video looks broken.
    pad = needle.length
    if needle.housing is not None:
        pad = max(pad, float(needle.housing.max()) * 0.6)
    lo = lo - pad * 0.5
    hi = hi + pad

    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")
    _equal_box(ax, lo, hi)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor("white")
        pane.pane.set_edgecolor(_GRID)
        pane._axinfo["grid"]["color"] = _GRID
    ax.set_xlabel("x (mm)"), ax.set_ylabel("y (mm)"), ax.set_zlabel("z (mm)")
    fmt = plt.FuncFormatter(lambda v, _: f"{v*1e3:.0f}")
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.zaxis.set_major_formatter(fmt)

    # Seed every collection with one degenerate segment: matplotlib 3.10's
    # add_collection3d autoscales from the segments and raises on an empty one.
    # They are cleared on the first draw.
    seed = [[path.xyz[0], path.xyz[0]]]
    wake = Line3DCollection(seed, colors=_WAKE, linewidths=2.6)
    ax.add_collection3d(wake)
    hot = Line3DCollection(seed, colors=_HOT, linewidths=3.4)
    ax.add_collection3d(hot)
    # The tool is drawn as a surface of revolution over its profile, so a luer
    # hub reads as a hub rather than a thick line. plot_surface cannot be
    # updated in place, so the artist is replaced each frame.
    tool_surf = [None]
    tip, = ax.plot([], [], [], color=_TIP, marker="o", ms=5.5, lw=0)
    box_lines = Line3DCollection(seed, colors=_NEEDLE, linewidths=1.1, alpha=0.55)
    ax.add_collection3d(box_lines)

    head = ax.set_title("", fontsize=11, loc="left")
    sub = fig.text(0.012, 0.022, "", fontsize=9, color="#6D6E71")

    # Wireframe edge list for the housing box (indices into Box.corners()).
    EDGES = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
             (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)]

    def draw(k):
        i = frames[k]
        segs = result.deposit.segments(upto=i)
        # Anything the tool came within `threshold` of gets drawn hot, so the
        # video marks the same rows the report counts.
        bad = {int(c) for c, d in zip(result.culprit, result.clearance)
               if c >= 0 and d < result.threshold}
        if len(segs):
            idx = np.arange(len(segs))
            mask = np.isin(idx, list(bad)) if bad else np.zeros(len(segs), bool)
            wake.set_segments(list(segs[~mask]))
            hot.set_segments(list(segs[mask]))
        else:
            wake.set_segments([])
            hot.set_segments([])

        rv = None if path.rotvec is None else path.rotvec[i]
        pose = needle.at(path.xyz[i], rv)
        if tool_surf[0] is not None:
            tool_surf[0].remove()
        X, Y, Z = needle.profile.surface(pose.tip, pose.R, n_theta=n_theta)
        tool_surf[0] = ax.plot_surface(X, Y, Z, color=_NEEDLE, alpha=0.95,
                                       linewidth=0, antialiased=True,
                                       shade=True, zorder=5)
        t = pose.tip
        tip.set_data([t[0]], [t[1]])
        tip.set_3d_properties([t[2]])

        if show_housing and pose.housing is not None:
            c = pose.housing.corners()
            box_lines.set_segments([[c[u], c[v]] for u, v in EDGES])
        else:
            box_lines.set_segments([])

        d = result.clearance[i]
        gap = "—" if not np.isfinite(d) else f"{d*1e3:.2f} mm"
        head.set_text(title or f"row {i+1}/{n}   ·   clearance {gap}")
        sub.set_text(f"{len(segs)} beads laid   ·   {result.status.upper()}   "
                     f"·   worst {result.worst*1e3:.2f} mm"
                     if np.isfinite(result.worst) else f"{len(segs)} beads laid")
        if spin:
            ax.view_init(elev=elev, azim=azim0 + spin * k / max(1, len(frames) - 1))
        if progress is not None:
            progress(k + 1, len(frames))
        return wake, hot, tip, box_lines

    ax.view_init(elev=elev, azim=azim0)
    anim = FuncAnimation(fig, draw, frames=len(frames), interval=1000 / fps,
                         blit=False)

    out = Path(out)
    writer = (PillowWriter(fps=fps) if out.suffix.lower() == ".gif"
              else FFMpegWriter(fps=fps, bitrate=2400))
    out.parent.mkdir(parents=True, exist_ok=True)
    anim.save(str(out), writer=writer)
    plt.close(fig)
    return out
