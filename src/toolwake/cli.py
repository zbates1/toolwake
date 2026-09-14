"""Command line: toolwake render <path> -o out.mp4"""
from __future__ import annotations

import argparse
import json
import sys

from .simulate import simulate
from .tool import Needle
from .toolpath import Toolpath


def _bar(label):
    def show(i, n):
        pct = 100.0 * i / n
        sys.stderr.write(f"\r  {label}: {pct:5.1f}%  ({i}/{n})")
        if i >= n:
            sys.stderr.write("\n")
        sys.stderr.flush()
    return show


def main(argv=None) -> int:
    p = argparse.ArgumentParser(
        prog="toolwake",
        description="Simulate a depositing tool and the wake it leaves behind.")
    p.add_argument("source", help="a .gcode file, or 'helix' for the built-in demo")
    p.add_argument("-o", "--out", default="wake.mp4",
                   help="video output (.mp4 needs ffmpeg, .gif does not)")
    p.add_argument("--inner-d", type=float, default=90e-6, help="needle ID, m")
    p.add_argument("--outer-d", type=float, default=None, help="needle OD, m")
    p.add_argument("--length", type=float, default=12.7e-3, help="needle length, m")
    p.add_argument("--housing", type=float, nargs=3, default=None,
                   metavar=("X", "Y", "Z"), help="end-effector extents, m")
    p.add_argument("--lag", type=int, default=8,
                   help="rows of immunity behind the nozzle")
    p.add_argument("--threshold", type=float, default=5e-4,
                   help="clearance below which a row is 'close', m")
    p.add_argument("--max-frames", type=int, default=600)
    p.add_argument("--fps", type=int, default=30)
    p.add_argument("--report", default=None, help="write the JSON report here")
    p.add_argument("--html", default=None,
                   help="write a standalone interactive viewer here")
    p.add_argument("--no-video", action="store_true", help="skip the video")
    a = p.parse_args(argv)

    path = (Toolpath.helix() if a.source == "helix"
            else Toolpath.from_gcode(a.source))
    print(f"  {path!r}", file=sys.stderr)

    needle = Needle(inner_d=a.inner_d, outer_d=a.outer_d, length=a.length,
                    housing=a.housing)
    res = simulate(path, needle, lag=a.lag, threshold=a.threshold,
                   progress=_bar("sweep"))
    rep = res.report()
    print(json.dumps(rep, indent=2))
    if a.report:
        with open(a.report, "w") as fh:
            json.dump(rep, fh, indent=2)

    if a.html:
        from .viewer import to_html
        print(f"  wrote {to_html(res, a.html)}", file=sys.stderr)

    if not a.no_video:
        from .render import animate
        out = animate(res, a.out, fps=a.fps, max_frames=a.max_frames,
                      progress=_bar("render"))
        print(f"  wrote {out}", file=sys.stderr)

    return 2 if rep["status"] == "collision" else 0


if __name__ == "__main__":
    raise SystemExit(main())
