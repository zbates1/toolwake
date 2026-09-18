# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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
