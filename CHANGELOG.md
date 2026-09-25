# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [0.3.0] — 2026-09-24

Performance. Verdicts are unchanged on every file tested; reported clearances
for DISTANT material are not, which is why this is a minor rather than a patch.

### Changed
- **`search` default 0.02 → 0.002 m.** Twenty millimetres of padding around
  every query meant that on a 10 mm part the query box covered the whole
  thing and the spatial hash could prune nothing: 949 candidate beads per
  query against a wake of 2 046, and 99.2% of the 30 million point-distance
  evaluations were against material nowhere near the tool.

  Detection is unaffected — 2 mm is still more than forty times the largest
  bead radius. What changes is that a clearance larger than `search` now
  reports `inf` instead of a number, so anything printing at a much larger
  scale wants this raised.
- **Grid cell size now follows the TOOL, not `search`.** They were tied
  (`cell = search / 4`), which welded two opposing costs to one knob:
  shrinking `search` cut candidates per query 36x (949 → 26) but grew the
  cells walked per query 1000x (8 → 7 913), because the cell shrank with it.
  Below about 5 mm the second swamped the first and asking for a tighter
  search radius made the whole sweep *slower*.

  Query boxes are a tool bounding box plus `search`, so the tool's own width
  is what the grid should be sized against. For a luer tip that is the 11 mm
  collar, giving 2.75 mm cells — within noise of the 2.5 mm that measured
  fastest. The housing is excluded: at 130 mm it would force cells so coarse
  that every needle query returned the whole wake.

Measured over the reference set, with every verdict and every worst clearance
identical to 0.2.3:

| toolpath | rows | 0.2.3 | 0.3.0 | |
|---|---|---|---|---|
| 360-180 SMOL Middle | 36 305 | 626 s | 252 s | 2.5x |
| 360-180 SMOL Middle (1) | 36 305 | 577 s | 261 s | 2.2x |
| test_infill | 2 674 | 7.3 s | 1.0 s | 7.2x |
| REDO WOW Slab | 2 674 | 7.8 s | 3.6 s | 2.2x |
| slab_infill | 2 674 | 7.2 s | 3.8 s | 1.9x |
| small_slab | 1 747 | 3.0 s | 2.3 s | 1.3x |
| LightPipeRobotTest | 1 170 | 0.8 s | 0.7 s | 1.1x |

Timings are from one laptop and vary with what else it is doing; treat the
ratios as indicative, not as a benchmark.

## [0.2.3] — 2026-09-18

Touching is not penetrating. Reported collision counts change.

### Fixed
- **Exact tangency was reported as a collision, at random.** With
  `bead_radius` and `bead_drop` both h/2 — the setting that makes a bead fill
  its layer, which 0.2.1 introduced — a same-layer neighbour passing under the
  nozzle wall sits at a clearance of EXACTLY zero. Which side of zero each one
  lands on is rounding, not geometry.

  On a real 36 305-row part that put **4 677 rows, 12.9% of the file**, on the
  wrong side and reported them as collisions. Every one was blamed on the
  cannula, every culprit bead was a same-layer neighbour 0.066–0.093 mm away
  laterally — under the wall, outside the bore — and every penetration depth
  printed as `-0.0000 mm`.

  `Result.n_hits` and `Result.status` now ignore anything shallower than
  `CONTACT_TOL` (1e-12 m), which is six orders of magnitude above
  double-precision noise on metre-scale coordinates and seven below anything
  this models. `report()` reports such a row as `0.0` rather than `-0.0`.
- **A zero threshold reported tangencies as "close".** With no warning band,
  a clearance of -1e-13 satisfied both "not penetrating" and "below the
  threshold", so a clean result came back as `warn`. A threshold of zero now
  means what it says: only penetration counts.

## [0.2.2] — 2026-09-18

Bounds the size of the viewer page. No reported number changes.

### Fixed
- **The viewer inlined every deposited bead, so the page grew without limit.**
  Frames were capped at 400; beads were not capped at all. A 2 674-row path
  already produces a 1.2 MB page, and a 100 000-row one would produce tens of
  megabytes — which is a problem the moment anyone raises their row budget and
  tries to embed the result in an iframe.

  `to_html_str` now takes `max_beads` (default 20 000) and decimates the drawn
  wake on the same rule as frames: sample evenly, and never drop a bead
  something collided with. Purely cosmetic — the report is computed from the
  full deposit, and nothing on the page claims a bead count.

## [0.2.1] — 2026-09-18

Finishes 0.2.0. That release fixed the tool's geometry; this one fixes the
material's. Reported clearances change again.

### Fixed
- **The bead was round at the bore diameter, so it was taller than its own
  layer.** A bead laid at layer height h is squashed to h and spreads
  sideways; kept round it pokes up through where the nozzle will sit on the
  next pass. On a 0.036 mm-layer slab with a 0.09 mm bore that was 498 rows
  reported, 391 of them blamed on material exactly one layer below.
- **The bead was centred on the toolpath, i.e. in the nozzle TIP's own plane.**
  Material leaving the bore fills the gap between the previous layer's top and
  the nozzle face, occupying [z - h, z] — its centre is h/2 down. Straddling
  the face's plane, any same-layer neighbour passing under the wall reported a
  collision of exactly one bead radius: tangent contacts dressed up as
  penetrations. 41 rows on the same slab, every one at exactly -bead_radius,
  which is what gave it away.

Rows reported as penetrating, across both releases:

| toolpath | 0.1.1 | 0.2.0 | 0.2.1 |
|---|---|---|---|
| LightPipeRobotTest | 7 | 0 | **0** |
| test_infill | 2 193 | 0 | **0** |
| small_slab | 1 224 | 498 | **0** |

Every G-code file in the reference set now reports clear.

### Added
- `simulate(bead_radius=..., bead_drop=...)` — the laid bead's size, and where
  it sits relative to the nozzle face. `bead_drop` is applied along the TOOL
  axis, not -Z, so it stays correct on a non-planar move. Explicit arguments
  rather than inferred from the path: only the caller knows whether a Z step
  is a layer or a ramp.

## [0.2.0] — 2026-09-18

Reported clearances change. The TOOL was modelled with two pieces of geometry
it does not have, and both made it collide with the print it was making. The
BEAD had two of its own; those are 0.2.1.

### Fixed
- **Tool sections were capsules, so the tool reached below its own tip.** A
  capsule's hemispherical cap bulges a full radius past its endpoint, so a
  0.095 mm-radius cannula occupied 0.095 mm of space beneath the nozzle — and
  the space beneath the nozzle is the part. On a 0.06 mm-layer print that
  swallowed the two layers underneath on nearly every row: 2 193 of 2 674 rows
  reported as collisions, worst -0.14 mm, every one blamed on the cannula
  (2 177 of them by material *below* the tip). Sections are now flat-ended
  `Cylinder`s. Beads stay capsules: extruded material does have rounded ends.
- **The cannula was modelled solid, so the bore counted as an obstacle.** A
  needle is a tube. Material under an open bore — the bead being laid, and the
  layer being laid onto — is in the hole, not against the wall. `Section` now
  takes `r_inner` and `Cylinder` an `inner_radius`. `blunt_cannula` already
  accepted `inner_d` and silently discarded it, which is how the bore went
  missing.

On the three reference toolpaths, rows reported as penetrating:

| toolpath | 0.1.1 | flat tip | + open bore |
|---|---|---|---|
| LightPipeRobotTest | 7 | 2 | **0** |
| test_infill | 2 193 | 16 | **0** |
| small_slab | 1 224 | 1 218 | 498 |

`small_slab`'s remaining 498 are the bead model, not the tool; fixed in 0.2.1.

### Added
- `Cylinder` — flat-ended, optionally hollow, exact signed distance.
- `toolwake --version` (also `-V` and `-v`). The CLI had no way to report its
  own version: `toolwake -v` failed with "the following arguments are
  required: source", which reads as a usage mistake rather than a missing
  feature.

## [0.1.1] — 2026-09-17

Performance only. Every report is byte-identical to 0.1.0 — verified on three
real toolpaths — because neither change touches what is computed, only how
much work is done to compute it.

### Fixed
- `Deposit._candidates` walked every grid cell in the query's bounding box.
  The tool's box includes the luer housing, 130 mm across, while the grid is
  sized for a 45 um bead, so each query swept roughly 5 000 mostly-empty cells:
  43.8 million dict lookups over one 1 747-row toolpath, and 86% of total
  runtime. It now enumerates whichever is smaller, the box's cells or the cells
  that actually hold material. Same set returned, so this is purely a cost
  choice.
- `Deposit` stored segments in Python lists and indexed them with
  `np.asarray(self._a)[cand]`, rebuilding an array of the entire wake on every
  query. Cost therefore tracked total deposit size rather than candidate
  count — quadratic over a sweep, and invisible behind a broad phase that was
  doing its job. Storage is now capacity-doubling arrays, so a query costs
  O(candidates). With candidates pinned at 2, a query went from 12.5 ms to
  0.57 ms at 16 000 segments.

Measured end to end, unchanged results:

| toolpath | rows | 0.1.0 | 0.1.1 |
|---|---|---|---|
| LightPipeRobotTest | 1 170 | 7.29 s | 0.42 s |
| small_slab | 1 747 | 16.80 s | 2.22 s |
| test_infill | 2 674 | 29.42 s | 6.72 s |

### Added
- `Deposit.last_probe` — cells inspected by the most recent query, so the cost
  invariant can be asserted without timing anything. Two regression tests use
  it; both were confirmed to fail against 0.1.0.

## [0.1.0] — 2026-09-14

First release.

### Added
- `Toolpath` — load from G-code (`M82`/`M83` aware, so print and travel moves
  are told apart correctly) or from `(N, 6)` pose files; `helix` and
  `conical_helix` fixtures.
- `Needle` / `ToolProfile` — tools as a coaxial stack of sections, the way CAM
  describes one. `Needle.luer()` models a luer-lock dispensing tip whose hub is
  roughly 50x the cannula diameter.
- `Deposit` — accumulating wake with a spatial hash for broad phase and exact
  capsule/box distance for narrow phase, plus a trailing-window rule so the
  tool never collides with the bead it is currently extruding.
- `simulate` — per-row clearance against everything laid earlier, naming which
  section of the tool was closest.
- `animate` — mp4/gif of the tool and its wake (needs matplotlib).
- `to_html` / `to_html_str` — a self-contained interactive viewer with no
  external requests; step between collisions, deep-link any frame.
- `ToolwakeSession` / `asgi_routes` — serve the viewer for embedding in another
  GUI, with a version poll so an iframe follows new results.
- `toolwake` CLI; exits 2 on a collision so it drops into a pre-flight gate.
