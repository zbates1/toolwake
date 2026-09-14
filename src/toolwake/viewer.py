"""A self-contained interactive HTML viewer for a simulation Result.

Writes ONE html file with no external requests — no CDN, no server, no Python
running. Open it, drag to rotate, and jump between collisions. That matters
because the people who most need to look at a failed clearance check are often
not the people with the environment installed.

Deliberately hand-rolled rather than plotly/three.js: the whole renderer is an
orthographic projection and a few hundred lines of canvas drawing, which keeps
the package dependency-free and the output openable offline forever.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .geometry import rotation_from_rotvec, rpy_from_rotation

__all__ = ["to_html", "to_html_str"]

_TEMPLATE = r"""<!doctype html>
<meta charset="utf-8">
<title>__TITLE__</title>
<style>
  :root {
    --ink:#1a1a1a; --muted:#6d6e71; --line:#d8dde4; --bg:#fff;
    --wake:#3fb8d4; --hot:#c41230; --tool:#2d2d8a; --tip:#ff9100;
    --ok:#1b7f4b; --warn:#c07a00;
  }
  * { box-sizing:border-box; }
  body { margin:0; font:14px/1.45 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;
         color:var(--ink); background:var(--bg); }
  .wrap { display:flex; height:100vh; }
  .stage { flex:1; position:relative; min-width:0; }
  canvas { display:block; width:100%; height:100%; cursor:grab; }
  canvas.drag { cursor:grabbing; }
  .side { width:290px; border-left:1px solid var(--line); padding:16px;
          overflow-y:auto; flex-shrink:0; }
  h1 { font-size:15px; margin:0 0 2px; }
  .sub { color:var(--muted); font-size:12px; margin-bottom:14px; }
  .lbl { color:var(--muted); font-size:11px; text-transform:uppercase;
         letter-spacing:.06em; margin:14px 0 4px; }
  .row { display:flex; gap:6px; margin:8px 0; }
  button { flex:1; padding:7px 4px; font:inherit; font-size:12px; cursor:pointer;
           background:#f2f4f7; border:1px solid var(--line); border-radius:7px; }
  button:hover:not(:disabled) { background:#e8ecf2; }
  button:disabled { opacity:.4; cursor:default; }
  button.hot { background:var(--hot); border-color:var(--hot); color:#fff; }
  input[type=range] { width:100%; }
  .hits { max-height:230px; overflow-y:auto; border:1px solid var(--line);
          border-radius:7px; margin-top:4px; }
  .hit { padding:5px 9px; font-size:12px; cursor:pointer; display:flex;
         justify-content:space-between; border-bottom:1px solid #eef1f5; }
  .hit:last-child{border-bottom:0}
  .hit:hover { background:#fdf3f5; }
  .hit.cur { background:var(--hot); color:#fff; }
  .hit .d { color:var(--hot); font-variant-numeric:tabular-nums; }
  .hit.cur .d { color:#fff; }
  .k { display:flex; justify-content:space-between; font-size:12px;
       padding:3px 0; border-bottom:1px dotted #eef1f5; }
  .k span:last-child { font-variant-numeric:tabular-nums; }
  .hud { position:absolute; left:14px; top:12px; font-size:13px; }
  .hud b { font-size:19px; }
  .hint { position:absolute; left:14px; bottom:10px; color:var(--muted);
          font-size:11px; }
</style>
<div class="wrap">
  <div class="stage">
    <canvas id="c"></canvas>
    <div class="hud"><b id="hrow"></b></div>
    <div class="hint">drag to rotate &middot; scroll to zoom</div>
  </div>
  <div class="side">
    <h1>__TITLE__</h1>
    <div class="sub">__SUBTITLE__</div>

    <div class="lbl">frame</div>
    <input type="range" id="scrub" min="0" value="0">
    <div class="row">
      <button id="prev" title="previous collision">&#9198; prev</button>
      <button id="play">&#9654; play</button>
      <button id="next" title="next collision">next &#9197;</button>
    </div>
    <div class="lbl">collisions (<span id="nhit"></span>)</div>
    <div class="hits" id="hits"></div>

    <div class="lbl">run</div>
    <div id="meta"></div>
  </div>
</div>
<script>
const D = __DATA__;

// ---- projection ---------------------------------------------------------
let yaw = -0.9, pitch = 0.42, zoom = 1, cx = 0, cy = 0;
const c = document.getElementById('c'), g = c.getContext('2d');
const lo = D.bounds[0], hi = D.bounds[1];
const mid = [0,1,2].map(i => (lo[i]+hi[i])/2);
const span = Math.max(hi[0]-lo[0], hi[1]-lo[1], hi[2]-lo[2]) || 1;

function proj(p) {
  const x = p[0]-mid[0], y = p[1]-mid[1], z = p[2]-mid[2];
  const cy_ = Math.cos(yaw), sy = Math.sin(yaw);
  const cp = Math.cos(pitch), sp = Math.sin(pitch);
  const X = x*cy_ - y*sy;
  const Y = x*sy + y*cy_;
  const s = Math.min(c.width, c.height) / (span*1.5) * zoom;
  return [c.width/2 + X*s + cx, c.height/2 - (z*cp - Y*sp)*s + cy];
}
function depth(p) {                       // painter's algorithm ordering
  const x = p[0]-mid[0], y = p[1]-mid[1], z = p[2]-mid[2];
  return (x*Math.sin(yaw) + y*Math.cos(yaw))*Math.cos(pitch) + z*Math.sin(pitch);
}

// ---- state --------------------------------------------------------------
let f = 0, playing = false, timer = null;
const N = D.frames.length;
const scrub = document.getElementById('scrub');
scrub.max = N - 1;

function resize() {
  const r = c.getBoundingClientRect(), dpr = window.devicePixelRatio || 1;
  c.width = r.width*dpr; c.height = r.height*dpr;
  g.setTransform(dpr,0,0,dpr,0,0);
  c.width = r.width*dpr; c.height = r.height*dpr;
  draw();
}

function draw() {
  const F = D.frames[f];
  g.setTransform(1,0,0,1,0,0);
  g.clearRect(0,0,c.width,c.height);

  // build plate
  g.strokeStyle = '#eef1f5'; g.lineWidth = 1;
  const r = span*0.75;
  for (let i=0;i<=8;i++){
    const t = -r + 2*r*i/8;
    let a = proj([mid[0]+t, mid[1]-r, lo[2]]), b = proj([mid[0]+t, mid[1]+r, lo[2]]);
    g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]); g.stroke();
    a = proj([mid[0]-r, mid[1]+t, lo[2]]); b = proj([mid[0]+r, mid[1]+t, lo[2]]);
    g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]); g.stroke();
  }

  // deposited beads laid up to this frame, back to front
  const vis = [];
  for (const s of D.beads) if (s[6] <= F.i) vis.push(s);
  vis.sort((p,q) => depth([p[0],p[1],p[2]]) - depth([q[0],q[1],q[2]]));
  g.lineCap = 'round';
  for (const s of vis) {
    const a = proj([s[0],s[1],s[2]]), b = proj([s[3],s[4],s[5]]);
    g.strokeStyle = s[7] ? getCss('--hot') : getCss('--wake');
    g.lineWidth = s[7] ? 3.4 : 2.4;
    g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]); g.stroke();
  }

  // the tool: one thick stroke per profile section, width from its radius
  const t = F.t, ax = F.a;
  const sc = Math.min(c.width, c.height) / (span*1.5) * zoom;
  let zoff = 0;
  for (const sec of D.profile) {
    const p0 = [t[0]+ax[0]*zoff, t[1]+ax[1]*zoff, t[2]+ax[2]*zoff];
    zoff += sec[0];
    const p1 = [t[0]+ax[0]*zoff, t[1]+ax[1]*zoff, t[2]+ax[2]*zoff];
    const a = proj(p0), b = proj(p1);
    g.strokeStyle = getCss('--tool');
    g.lineWidth = Math.max(1.5, 2*Math.max(sec[1],sec[2])*sc);
    g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]); g.stroke();
  }
  const tp = proj(t);
  g.fillStyle = getCss('--tip');
  g.beginPath(); g.arc(tp[0],tp[1], 4.5, 0, 6.2832); g.fill();

  // housing wireframe
  if (F.box) {
    g.strokeStyle = getCss('--tool'); g.globalAlpha = .45; g.lineWidth = 1.2;
    for (const e of D.edges) {
      const a = proj(F.box[e[0]]), b = proj(F.box[e[1]]);
      g.beginPath(); g.moveTo(a[0],a[1]); g.lineTo(b[0],b[1]); g.stroke();
    }
    g.globalAlpha = 1;
  }

  document.getElementById('hrow').textContent = `row ${F.i+1} / ${D.rows}`;
  scrub.value = f;
  markList();
}
function getCss(v){ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }

// ---- collision navigation ----------------------------------------------
const hits = D.collisions;
function nearestHit(dir) {
  if (!hits.length) return null;
  if (dir > 0) { for (const h of hits) if (h > f) return h; return hits[0]; }
  for (let i=hits.length-1;i>=0;i--) if (hits[i] < f) return hits[i];
  return hits[hits.length-1];
}
function go(i, quiet){
  f = Math.max(0, Math.min(N-1, i));
  // Reflect the frame in the URL so a specific collision can be sent to
  // someone else: "open this and look at row 506" becomes a link.
  if (!quiet) history.replaceState(null, '', '#f=' + D.frames[f].i);
  draw();
}

document.getElementById('next').onclick = () => { const h = nearestHit(1); if (h!==null) go(h); };
document.getElementById('prev').onclick = () => { const h = nearestHit(-1); if (h!==null) go(h); };
document.getElementById('play').onclick = function () {
  playing = !playing;
  this.innerHTML = playing ? '&#10073;&#10073; pause' : '&#9654; play';
  if (playing) timer = setInterval(() => go(f >= N-1 ? 0 : f+1), 1000/D.fps);
  else clearInterval(timer);
};
scrub.oninput = e => go(+e.target.value);
document.onkeydown = e => {
  if (e.key === 'ArrowRight') go(f+1);
  if (e.key === 'ArrowLeft') go(f-1);
  if (e.key === 'n') document.getElementById('next').click();
  if (e.key === 'p') document.getElementById('prev').click();
};

// ---- collision list -----------------------------------------------------
const hitsEl = document.getElementById('hits');
document.getElementById('nhit').textContent = hits.length;
if (!hits.length) hitsEl.innerHTML =
  '<div class="hit" style="cursor:default;color:var(--muted)">none — path is clear</div>';
for (const h of hits) {
  const d = document.createElement('div');
  d.className = 'hit'; d.dataset.f = h;
  d.innerHTML = `<span>row ${D.frames[h].i+1}</span>`
              + `<span class="d">${(D.frames[h].c*1000).toFixed(2)} mm</span>`;
  d.onclick = () => go(h);
  hitsEl.appendChild(d);
}
function markList(){
  for (const el of hitsEl.children)
    el.classList.toggle('cur', +el.dataset.f === f);
}

const meta = document.getElementById('meta');
for (const [k,v] of D.meta)
  meta.insertAdjacentHTML('beforeend', `<div class="k"><span>${k}</span><span>${v}</span></div>`);

// ---- mouse --------------------------------------------------------------
let drag = null;
c.onmousedown = e => { drag = [e.clientX, e.clientY]; c.classList.add('drag'); };
window.onmouseup = () => { drag = null; c.classList.remove('drag'); };
window.onmousemove = e => {
  if (!drag) return;
  yaw += (e.clientX-drag[0])*0.01;
  pitch = Math.max(-1.5, Math.min(1.5, pitch + (e.clientY-drag[1])*0.01));
  drag = [e.clientX, e.clientY]; draw();
};
c.onwheel = e => { e.preventDefault();
  zoom = Math.max(0.25, Math.min(8, zoom * (e.deltaY < 0 ? 1.1 : 1/1.1))); draw(); };

// Deep link: #worst, or #f=<row index as printed in the report>.
(function () {
  const h = decodeURIComponent(location.hash || '').replace('#', '');
  if (h === 'worst') { f = D.worst_frame; return; }
  const m = h.match(/^f=(\d+)$/);
  if (!m) return;
  const want = +m[1];
  let best = 0, bd = Infinity;            // nearest stored frame to that row
  D.frames.forEach((fr, k) => { const d = Math.abs(fr.i - want);
                                if (d < bd) { bd = d; best = k; } });
  f = best;
})();

window.onresize = resize;
resize();
</script>
"""

def _meta_rows(result, n, rep):
    """The summary table: how big the cloud is, how far the tool turns, what hit.

    Roll/pitch/yaw come from the per-row orientation in the ZYX convention a
    teach pendant reports, so the range here is directly comparable to what an
    operator reads off the machine. A planar path has no orientation, so the
    ranges are reported as flat zero rather than omitted — "0 to 0" is a fact
    about the path, not missing data.
    """
    path = result.path
    if path.rotvec is None:
        rpy = np.zeros((len(path), 3))
    else:
        rpy = np.array([rpy_from_rotation(rotation_from_rotvec(rv))
                        for rv in path.rotvec])
    deg = np.degrees(rpy)
    rows = [("point cloud size", f"{n:,}")]
    for k, name in enumerate(("roll", "pitch", "yaw")):
        rows.append((f"{name} min / max",
                     f"{deg[:, k].min():.1f}° / {deg[:, k].max():.1f}°"))
    rows.append(("colliding frames", f"{rep['rows_penetrating']:,}"))
    return rows


_EDGES = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
          (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)]


def to_html(result, out="wake.html", **kw):
    """Write the viewer to a file. See `to_html_str` for the options."""
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(to_html_str(result, **kw), encoding="utf-8")
    return out


def to_html_str(result, *, max_frames=400, fps=20, title=None,
                show_housing=True, live_poll=None):
    """Write a standalone interactive viewer for `result`.

    Buttons step between collision frames, the list jumps to any of them, and
    "jump to worst" goes straight to the deepest penetration — which is usually
    the only frame anyone actually wants to see.

    The current frame is mirrored into the URL hash, so `file.html#f=506` or
    `file.html#worst` opens directly on a given row. That turns "look at the
    collision around row 506" into something sendable.

    Args:
        result: a `Result` from `toolwake.simulate`.
        out: destination `.html`.
        max_frames: cap on stored frames; the path is subsampled to fit. Every
            collision frame is kept regardless, so nothing navigable is lost.
        fps: playback rate for the play button.
        show_housing: include the end-effector wireframe.

    Args (continued):
        live_poll: when set, a URL the page polls for a version number,
            reloading itself when it changes. This is what lets an embedded
            iframe follow new results without the host page doing anything.

    Returns:
        The HTML, as a string.
    """
    path, needle, dep = result.path, result.needle, result.deposit
    n = len(path)

    # Subsample for size, but never drop a collision: the whole point of the
    # viewer is getting to those frames.
    step = max(1, int(np.ceil(n / max_frames)))
    keep = set(range(0, n, step)) | {n - 1}
    keep |= {int(i) for i in np.nonzero(result.clearance < 0)[0]}
    rows = sorted(keep)
    index_of = {r: k for k, r in enumerate(rows)}

    hot = {int(c) for c, d in zip(result.culprit, result.clearance)
           if c >= 0 and d < result.threshold}
    segs = dep.segments()
    t_of = np.asarray(dep._t) if len(dep) else np.zeros(0, dtype=int)
    beads = [[*map(float, s[0]), *map(float, s[1]), int(t_of[j]), 1 if j in hot else 0]
             for j, s in enumerate(segs)]

    frames = []
    for r in rows:
        rv = None if path.rotvec is None else path.rotvec[r]
        pose = needle.at(path.xyz[r], rv)
        d = result.clearance[r]
        frames.append({
            "i": int(r),
            "t": [float(v) for v in pose.tip],
            "a": [float(v) for v in pose.axis],
            "c": None if not np.isfinite(d) else float(d),
            "src": str(result.source[r]),
            "box": ([[float(v) for v in p] for p in pose.housing.corners()]
                    if (show_housing and pose.housing is not None) else None),
        })

    lo, hi = path.bounds()
    if len(segs):
        lo = np.minimum(lo, segs.reshape(-1, 3).min(axis=0))
        hi = np.maximum(hi, segs.reshape(-1, 3).max(axis=0))
    hi = hi + needle.length * 0.6

    collisions = [index_of[int(i)] for i in np.nonzero(result.clearance < 0)[0]
                  if int(i) in index_of]
    worst_row = int(np.argmin(np.where(np.isfinite(result.clearance),
                                       result.clearance, np.inf)))
    worst_frame = index_of.get(worst_row, 0)

    rep = result.report()
    data = {
        "rows": n,
        "fps": fps,
        "threshold": float(result.threshold),
        "bounds": [[float(v) for v in lo], [float(v) for v in hi]],
        "profile": [[float(s.length), float(s.r0), float(s.r_top)]
                    for s in needle.profile.sections],
        "beads": beads,
        "frames": frames,
        "collisions": collisions,
        "worst_frame": worst_frame,
        "edges": _EDGES,
        "meta": _meta_rows(result, n, rep),
    }

    html = (_TEMPLATE
            .replace("__DATA__", json.dumps(data, separators=(",", ":")))
            .replace("__TITLE__", title or "Robot arm printing needle collision checker")
            .replace("__SUBTITLE__",
                     f"Point cloud size {n:,}  ·  "
                     f"{'6DoF' if not path.is_planar else 'planar'}"))

    if live_poll and "</script>" in html:
        html = html.replace(
            "</script>",
            _LIVE_JS.replace("__URL__", live_poll) + "</script>", 1)
    return html


# Appended only when served: poll a version endpoint and reload on change, so
# an embedded iframe follows new checks without the host page doing anything.
_LIVE_JS = """
(function () {
  let seen = null;
  setInterval(async () => {
    try {
      const r = await fetch('__URL__', {cache: 'no-store'});
      const v = (await r.json()).version;
      if (seen === null) { seen = v; return; }
      if (v !== seen) location.reload();
    } catch (e) { /* host gone; keep showing the last result */ }
  }, 1500);
})();
"""
