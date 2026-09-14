# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[Semantic Versioning](https://semver.org/).

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
