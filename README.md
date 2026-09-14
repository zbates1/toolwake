# toolwake

Simulate a depositing tool travelling a toolpath, accumulate the material it
leaves in its wake, and check the tool against that wake.

Built for non-planar bioprinting, where a tall part and a wide end effector can
collide with material the printer laid down minutes earlier — a failure that
ordinary slicer previews and robot self-collision checks both miss.

The video and the numbers come from the same model, so the picture cannot
disagree with the report.

```python
from toolwake import Toolpath, Needle, simulate, animate

path   = Toolpath.helix(radius=0.018, pitch=0.0035, turns=7)
needle = Needle(inner_d=90e-6, outer_d=0.4e-3, length=12.7e-3,
                housing=(0.130, 0.075, 0.072))

res = simulate(path, needle, lag=10)
print(res.report())
# {'status': 'ok', 'rows': 770, 'min_clearance_mm': 3.2626, ...}

animate(res, "wake.mp4")        # video
to_html(res, "wake.html")       # interactive, self-contained
```

The HTML viewer is one file with no external requests — drag to rotate, scrub
the timeline, and step between collisions. `wake.html#worst` opens on the
deepest penetration, and every frame is linkable as `#f=<row>`, so "look at the
collision around row 506" becomes something you can send.

Or from the shell:

```bash
toolwake helix -o wake.mp4
toolwake part.gcode --housing 0.130 0.075 0.072 --report clearance.json
toolwake part.gcode --html check.html --no-video
```

The CLI exits `2` on a collision, so it drops straight into a pre-flight gate.

## Install

```bash
pip install toolwake            # core: numpy only
pip install toolwake[render]    # + matplotlib for video
                                # the HTML viewer needs neither
```

`.mp4` output needs ffmpeg on PATH. `.gif` does not.

## Why it is built this way

**The tool is a profile, the way CAM describes one.** A dispensing needle is
not a single cylinder: a luer-lock tip is a thin cannula stepping up to a hub
roughly 50x its diameter, then a locking collar. Sections stack from the tip
upward, each a cylinder or a taper, so the model has the same holder / shank /
flute structure a CAM simulation draws.

```python
Needle.luer(inner_d=90e-6)
# ToolProfile(4 sections: cannula, taper, hub, luer collar;
#             len=32.7 mm, r_max=5.50 mm)
```

Reports name the section that was closest, because "hub" and "cannula" call for
different fixes.

**Primitives, not meshes.** Each section is a capsule and the end effector is an
oriented box, so the distance functions are exact and there is no mesh
dependency. Tapers become capsules at their LARGEST radius — deliberately
conservative, since a false stop costs a reprint and a missed collision costs
the part. Pass `slices=` to follow a taper more tightly.

**Model the housing, not just the tip.** The nozzle crossing an earlier bead is
the failure people picture and the rarer one. On a tall part the housing is down
inside the structure while the tip is somewhere harmless. `test_housing_catches_
what_the_needle_misses` is that case in one test: a thin needle reports clear,
a 130 mm body reports a collision, on the same path.

**Non-planar is orientation, not just Z.** A 6-DOF path carries a rotation
vector per row, read the way a UR reports pose, so rows move between the robot
and here unchanged. `conical_helix` leans the tool to the wall's own slope;
`with_tool_axis` turns any set of surface normals into orientations.

Orientation only changes the answer when the closest approach is somewhere
other than the tip — lean the tool and the tip has not moved, but the hub has
swung several millimetres. That is the whole reason the body is modelled.

**Deposit after measuring.** Each row is checked against the wake *before* its
own material is added, and `lag` extends that immunity a few rows back —
otherwise the tool always collides with the bead it is extruding right now.

**Only print moves leave a wake.** Travel moves deposit nothing. Treating them
as material fills the part with phantom obstacles and every later move reports a
false collision. `Toolpath.from_gcode` handles both `M82` and `M83`, because
whether a row extrudes depends on it and slicers disagree — Cura emits `M82`,
PrusaSlicer `M83`.

**Two-phase queries.** A uniform spatial hash narrows the candidate beads, then
exact capsule/box distance measures against those only. A plain occupancy grid
would give a binary answer quantised to the voxel size; this reports a real
distance. The grid is sized to the *query* radius rather than the bead — getting
that wrong made an early version 16x slower for identical results.

## Embedding in another GUI

`ToolwakeSession` holds the latest result and hands it out as HTML, so a host
application can iframe it and let it follow new checks on its own.

```python
from toolwake import ToolwakeSession, Needle

session = ToolwakeSession(needle=Needle.luer(inner_d=90e-6), lag=10)
url = session.serve(port=7100)        # -> http://127.0.0.1:7100
...
session.run(toolpath)                 # the iframe reloads itself
```

The page polls a version endpoint, so re-checking a toolpath updates the view
without the host telling it anything.

If the host already runs a web framework, skip the extra port and mount the
routes — same bytes, no second socket, no cross-origin handling:

```python
from fastapi import Response
from toolwake.serve import asgi_routes

for path, fn in asgi_routes(session):
    def make(fn=fn):
        def endpoint():
            body, ctype = fn()
            return Response(body, media_type=ctype)
        return endpoint
    app.get("/api/toolwake" + path)(make())
```

Then point an iframe at `/api/toolwake/view`. Routes are `/view`, `/version`
and `/report.json`.

> On some Windows machines a local security product stalls large responses over
> loopback — measured here as 1 KB instant and 83 KB failing under a bare
> stdlib handler with no toolwake code in the path. If `serve()` hangs, mount
> the routes instead.

## API

| | |
|---|---|
| `Toolpath.from_gcode(path)` | parse G0/G1, infer print vs travel from E |
| `Toolpath.from_arrays(xyz, kinds, rotvec)` | bring your own, 3-DOF or 6-DOF |
| `Toolpath.from_poses(path)` | (N, 6) `x y z rx ry rz` pose file |
| `Toolpath.with_tool_axis(axes)` | set orientation from per-row tool directions |
| `Toolpath.helix(...)` | vertical-tool helix fixture |
| `Toolpath.conical_helix(...)` | flaring wall; tool leans to stay normal |
| `Needle(inner_d, outer_d, length, housing)` | plain cannula |
| `Needle.luer(inner_d, ...)` | full luer tip: cannula, taper, hub, collar |
| `ToolProfile([Section(...), ...])` | build any coaxial tool |
| `simulate(path, needle, lag, threshold)` | → `Result` |
| `Result.report()` | JSON-ready clearance summary |
| `animate(result, "out.mp4")` | video of the tool and its wake |
| `to_html(result, "out.html")` | interactive viewer, jumps between collisions |
| `ToolwakeSession(...)` | re-checkable session, serves HTML for embedding |
| `asgi_routes(session)` | route table to mount into an existing app |
| `Deposit`, `Capsule`, `Box` | the pieces, usable directly |

Distances are **metres** throughout; the report converts to mm.

`Result.status` is `ok` / `warn` / `collision`, matching the usual clearance
banding: penetration is a collision, anything inside `threshold` is a warning.

## Limits

- **Per-row sampling, not swept volumes.** A thin tool moving fast past a thin
  wall can pass between sampled rows. Densify the path or add continuous
  collision detection if that matters for your geometry.
- **Beads are straight segments** between consecutive rows, with a constant
  radius. Die swell, sag and variable extrusion width are not modelled.
- **No process physics** — nothing here knows about cure, flow or adhesion.
- **Preset hub dimensions are nominal.** `luer_taper_tip` defaults approximate a
  34G half-inch tip; hubs vary between manufacturers. Measure the tip you run
  and pass the real numbers.
- `search` bounds how far a query looks. Beads beyond it are not measured and
  the clearance returns `inf`. Raise it for a large housing.

## Development

```bash
pip install -e .[dev]
pytest
```

## Licence

MIT
