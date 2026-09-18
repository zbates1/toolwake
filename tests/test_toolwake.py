"""Tests for toolwake.

The important ones are not "does it run" but "does it catch a collision it
should, and stay quiet on one it shouldn't". Both directions are checked,
because a clearance test that never fires is indistinguishable from one that
works until the day it matters.
"""
import time
import numpy as np
import pytest

from toolwake import (Box, Capsule, Deposit, Needle, Toolpath, simulate,
                      rotation_from_rotvec)
from toolwake.geometry import rotvec_from_axis
from toolwake.toolpath import PRINT, TRAVEL


# ------------------------------------------------------------------ geometry
def test_capsule_distance_is_exact():
    c = Capsule([0, 0, 0], [0, 0, 1], radius=0.1)
    assert c.distance([[0.5, 0, 0.5]])[0] == pytest.approx(0.4)   # radial
    assert c.distance([[0, 0, 0.5]])[0] == pytest.approx(-0.1)    # on the axis
    # Beyond the end cap the distance is measured to the cap, not the line.
    assert c.distance([[0, 0, 2.0]])[0] == pytest.approx(0.9)


def test_capsule_uses_segment_not_infinite_line():
    c = Capsule([0, 0, 0], [0, 0, 1], radius=0.0)
    far = c.distance([[0, 0, 11.0]])[0]
    assert far == pytest.approx(10.0), "distance must be to the end, not the axis"


def test_box_distance_inside_and_out():
    b = Box([0, 0, 0], [1, 1, 1])
    assert b.distance([[2, 0, 0]])[0] == pytest.approx(1.0)
    assert b.distance([[2, 2, 0]])[0] == pytest.approx(np.sqrt(2))   # corner
    assert b.distance([[0, 0, 0]])[0] == pytest.approx(-1.0)         # centre


def test_box_rotation_is_applied():
    R = rotation_from_rotvec([0, 0, np.pi / 4])
    b = Box([0, 0, 0], [1, 0.1, 1], R)
    # The long axis now points along the diagonal, so a point out along it is
    # closer than the same distance along x would be for an unrotated box.
    on_axis = b.distance([[1.0 / np.sqrt(2), 1.0 / np.sqrt(2), 0]])[0]
    assert on_axis == pytest.approx(0.0, abs=1e-9)


def test_rotvec_identity_and_known_angle():
    assert np.allclose(rotation_from_rotvec([0, 0, 0]), np.eye(3))
    R = rotation_from_rotvec([0, 0, np.pi / 2])
    assert np.allclose(R @ [1, 0, 0], [0, 1, 0], atol=1e-12)


# ------------------------------------------------------------------- deposit
def test_deposit_respects_time_order():
    d = Deposit(bead_radius=0.0)
    d.add([0, 0, 0], [1, 0, 0], frame=10)
    probe = Capsule([0.5, 0, 0], [0.5, 0, 1], radius=0.0)
    # Laid at 10: invisible to an earlier frame, visible to a later one.
    assert d.clearance(probe, frame=5, lag=0)[0] == float("inf")
    assert d.clearance(probe, frame=50, lag=0)[0] == pytest.approx(0.0)


def test_deposit_lag_grants_immunity():
    d = Deposit(bead_radius=0.0)
    d.add([0, 0, 0], [1, 0, 0], frame=10)
    probe = Capsule([0.5, 0, 0], [0.5, 0, 1], radius=0.0)
    assert d.clearance(probe, frame=12, lag=8)[0] == float("inf")
    assert d.clearance(probe, frame=30, lag=8)[0] == pytest.approx(0.0)


def test_deposit_subtracts_bead_radius():
    d = Deposit(bead_radius=0.1)
    d.add([0, 0, 0], [1, 0, 0], frame=0)
    probe = Capsule([0.5, 0, 1.0], [0.5, 0, 2.0], radius=0.0)
    # search must exceed the gap: these are unit-scale shapes, well outside the
    # 2 cm default tuned for a bioprinter. Raising it here is the point — the
    # default silently returns inf, which test_deposit_far_away_is_out_of_range
    # pins down deliberately.
    assert d.clearance(probe, frame=99, lag=0, search=2.0)[0] == pytest.approx(0.9)


def test_deposit_far_away_is_out_of_range():
    d = Deposit(bead_radius=0.0)
    d.add([0, 0, 0], [1, 0, 0], frame=0)
    probe = Capsule([0, 0, 50.0], [0, 0, 51.0], radius=0.0)
    assert d.clearance(probe, frame=99, lag=0, search=0.02)[0] == float("inf")


# ------------------------------------------------------------------ toolpath
def test_only_print_rows_leave_a_wake():
    xyz = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]], float)
    path = Toolpath(xyz, kinds=[PRINT, PRINT, TRAVEL, PRINT])
    res = simulate(path, Needle(inner_d=1e-3, length=1e-3), lag=0)
    # 3 moves follow row 0; one is travel, so only 2 beads should exist.
    assert len(res.deposit) == 2


def test_gcode_roundtrip(tmp_path):
    g = tmp_path / "t.gcode"
    g.write_text("\n".join([
        "M82", "G1 X0 Y0 Z0 F600",
        "G1 X10 Y0 Z0 E1",      # extrude -> print
        "G1 X20 Y0 Z0",         # no E -> travel
        "G1 X30 Y0 Z0 E2",      # absolute E rises -> print
    ]))
    p = Toolpath.from_gcode(g)
    assert len(p) == 3
    assert list(p.kinds) == [PRINT, TRAVEL, PRINT]
    assert p.xyz[0] == pytest.approx([0.010, 0, 0])     # mm -> m


def test_gcode_relative_extrusion_m83(tmp_path):
    g = tmp_path / "rel.gcode"
    g.write_text("\n".join([
        "M83", "G1 X0 Y0 Z0 F600",
        "G1 X10 Y0 Z0 E0.5",
        "G1 X20 Y0 Z0 E0.5",
    ]))
    p = Toolpath.from_gcode(g)
    assert list(p.kinds) == [PRINT, PRINT], "M83 deltas are extrusions, not totals"


# ------------------------------------------------------------------ simulate
def test_clean_path_reports_ok():
    """A single straight line cannot collide with its own wake."""
    xyz = np.column_stack([np.linspace(0, 0.05, 60), np.zeros(60), np.zeros(60)])
    res = simulate(Toolpath(xyz), Needle(inner_d=90e-6), lag=4)
    assert res.status == "ok"
    assert res.n_hits == 0


def test_revisiting_a_layer_is_caught():
    """Come back along a line already printed, at the same height.

    The needle has real thickness, so re-treading its own bead must register.
    """
    fwd = np.column_stack([np.linspace(0, 0.05, 40), np.zeros(40), np.zeros(40)])
    back = fwd[::-1].copy()
    path = Toolpath(np.vstack([fwd, back]))
    res = simulate(path, Needle(inner_d=1e-3, outer_d=2e-3, length=5e-3), lag=4)
    assert res.status == "collision", res.report()
    assert res.n_hits > 0


def test_housing_catches_what_the_needle_misses():
    """The whole argument for modelling the body, in one test.

    A tall wall is printed, then the tool moves alongside it, far enough away
    that the thin needle is clear but close enough that a wide housing is not.
    """
    z = np.linspace(0, 0.06, 120)
    wall = np.column_stack([np.zeros(120), np.zeros(120), z])
    # Travel back down 25 mm to the side: clear for a 1 mm needle, not for a
    # 130 mm-wide body.
    aside = np.column_stack([np.full(40, 0.025), np.zeros(40),
                             np.linspace(0.06, 0.0, 40)])
    path = Toolpath(np.vstack([wall, aside]),
                    kinds=[PRINT] * 120 + [TRAVEL] * 40)

    thin = simulate(path, Needle(inner_d=1e-3, outer_d=2e-3, length=6e-3), lag=4)
    wide = simulate(path, Needle(inner_d=1e-3, outer_d=2e-3, length=6e-3,
                                 housing=(0.130, 0.075, 0.072)), lag=4)

    assert thin.worst > 0, "a thin needle should stay clear here"
    assert wide.n_hits > 0, "the housing should collide where the needle does not"
    assert wide.worst < thin.worst


def test_report_is_json_shaped():
    res = simulate(Toolpath.helix(turns=1.0, per_turn=40), Needle())
    rep = res.report()
    for k in ("status", "rows", "deposited_segments", "min_clearance_mm",
              "rows_penetrating", "threshold_mm", "lag_rows"):
        assert k in rep
    import json
    json.loads(json.dumps(rep))          # must survive a round trip


def test_helix_fixture_is_sane():
    p = Toolpath.helix(radius=0.02, pitch=0.004, turns=3, per_turn=60)
    assert len(p) == 180
    lo, hi = p.bounds()
    assert hi[2] - lo[2] == pytest.approx(0.012, rel=1e-3)   # 3 turns x 4 mm


# ------------------------------------------------------------------- profile
def test_profile_stacks_tip_first():
    from toolwake import Section, ToolProfile
    p = ToolProfile([Section(0.010, 0.001, name="a"),
                     Section(0.020, 0.004, name="b")])
    assert p.length == pytest.approx(0.030)
    assert p.tip_radius == pytest.approx(0.001)    # tip radius is the FIRST
    assert p.r_max == pytest.approx(0.004)


def test_profile_capsules_are_placed_along_the_axis():
    from toolwake import Section, ToolProfile
    p = ToolProfile([Section(0.010, 0.001, name="a"),
                     Section(0.020, 0.004, name="b")])
    caps = p.capsules([0, 0, 0], [0, 0, 1])
    assert [n for n, _ in caps] == ["a", "b"]
    assert caps[0][1].a == pytest.approx([0, 0, 0])
    assert caps[0][1].b == pytest.approx([0, 0, 0.010])   # b starts where a ends
    assert caps[1][1].a == pytest.approx([0, 0, 0.010])
    assert caps[1][1].b == pytest.approx([0, 0, 0.030])


def test_taper_capsule_is_conservative():
    """A cone becomes a capsule at its LARGEST radius, never its smallest."""
    from toolwake import Section, ToolProfile
    p = ToolProfile([Section(0.010, 0.001, 0.005, name="cone")])
    (_, cap), = p.capsules([0, 0, 0], [0, 0, 1], slices=1)
    assert cap.radius == pytest.approx(0.005)


def test_taper_slicing_tightens_the_fit():
    from toolwake import Section, ToolProfile
    p = ToolProfile([Section(0.010, 0.001, 0.005, name="cone")])
    coarse = p.capsules([0, 0, 0], [0, 0, 1], slices=1)
    fine = p.capsules([0, 0, 0], [0, 0, 1], slices=5)
    assert len(fine) == 5
    # The slice nearest the tip must be thinner than the single fat capsule.
    assert fine[0][1].radius < coarse[0][1].radius


def test_luer_hub_dwarfs_the_cannula():
    """The reason the profile exists, as a number."""
    from toolwake import luer_taper_tip
    p = luer_taper_tip()
    assert p.r_max / p.tip_radius > 20


def test_luer_hub_collides_where_cannula_does_not():
    """Same path, same tip position — only the hub reaches the wall.

    A bare cannula passes; the luer tip does not. That gap is exactly the
    failure a tip-only check waves through.
    """
    n = 120
    wall = np.column_stack([np.zeros(n), np.zeros(n), np.linspace(0, 0.05, n)])
    # Come back down 3 mm to the side: clear for a 0.19 mm cannula, not for a
    # 9 mm hub.
    aside = np.column_stack([np.full(50, 0.003), np.zeros(50),
                             np.linspace(0.05, 0.02, 50)])
    path = Toolpath(np.vstack([wall, aside]),
                    kinds=[PRINT] * n + [TRAVEL] * 50)

    bare = simulate(path, Needle(inner_d=90e-6, outer_d=0.19e-3,
                                 length=12.7e-3), lag=4)
    luer = simulate(path, Needle.luer(inner_d=90e-6), lag=4)

    assert bare.n_hits == 0, "a bare cannula should clear a 3 mm gap"
    assert luer.n_hits > 0, "a 9 mm hub should not"


def test_report_names_the_offending_part():
    """'hub' and 'cannula' need different fixes, so the report must say which."""
    n = 120
    wall = np.column_stack([np.zeros(n), np.zeros(n), np.linspace(0, 0.05, n)])
    aside = np.column_stack([np.full(50, 0.003), np.zeros(50),
                             np.linspace(0.05, 0.02, 50)])
    path = Toolpath(np.vstack([wall, aside]), kinds=[PRINT] * n + [TRAVEL] * 50)
    res = simulate(path, Needle.luer(inner_d=90e-6), lag=4)
    assert res.report()["min_clearance_source"] in {"hub", "taper", "luer collar"}


# ----------------------------------------------------------------- nonplanar
def test_rotvec_from_axis_round_trips():
    for u in ([0, 0, 1], [1, 0, 0], [0, 1, 0], [1, 1, 1], [-0.3, 0.7, 0.2]):
        got = rotation_from_rotvec(rotvec_from_axis(u)) @ [0, 0, 1]
        want = np.asarray(u, float) / np.linalg.norm(u)
        assert got == pytest.approx(want, abs=1e-9)


def test_rotvec_from_axis_handles_antiparallel():
    """The 180 deg case has no unique axis and must not divide by zero."""
    got = rotation_from_rotvec(rotvec_from_axis([0, 0, -1])) @ [0, 0, 1]
    assert got == pytest.approx([0, 0, -1], abs=1e-9)


def test_with_tool_axis_sets_orientation():
    xyz = np.zeros((4, 3))
    axes = np.tile([1.0, 0.0, 0.0], (4, 1))
    p = Toolpath(xyz).with_tool_axis(axes)
    assert not p.is_planar
    for i in range(4):
        got = rotation_from_rotvec(p.rotvec[i]) @ [0, 0, 1]
        assert got == pytest.approx([1, 0, 0], abs=1e-9)


def test_needle_follows_the_tool_axis():
    """A tilted row must place the tool body sideways, not just record a number."""
    n = Needle(inner_d=1e-3, outer_d=2e-3, length=0.010)
    pose = n.at([0, 0, 0], rotvec=rotvec_from_axis([1, 0, 0]))
    assert pose.needle.b == pytest.approx([0.010, 0, 0], abs=1e-9)


def test_conical_helix_lean_matches_the_wall():
    p = Toolpath.conical_helix(r0=0.008, r1=0.018, height=0.030,
                               turns=4, per_turn=60)
    assert not p.is_planar
    expect = np.arctan((0.018 - 0.008) / 0.030)
    for i in (0, len(p) // 2, len(p) - 1):
        ax = rotation_from_rotvec(p.rotvec[i]) @ [0, 0, 1]
        lean = np.arccos(np.clip(ax[2], -1, 1))
        assert lean == pytest.approx(expect, abs=1e-6)


def test_conical_helix_tilt_zero_is_vertical():
    p = Toolpath.conical_helix(turns=2, per_turn=40, tilt=0.0)
    for i in (0, len(p) - 1):
        ax = rotation_from_rotvec(p.rotvec[i]) @ [0, 0, 1]
        assert ax == pytest.approx([0, 0, 1], abs=1e-9)


def test_orientation_changes_clearance_of_the_body():
    """The point of non-planar support, as a measurement.

    Identical tip position, two tool leans, against a wall beside it. Any
    difference can only come from the BODY swinging — which is exactly what
    orientation is tracked for.

    Note what this does NOT claim: if the closest approach happens at the tip,
    leaning changes nothing, because the tip has not moved. An earlier version
    of this test asserted otherwise on a conical helix and failed, both runs
    returning 4.464 mm — correctly, because the cannula tip was binding.
    """
    dep = Deposit(bead_radius=1e-4)
    for k, z in enumerate(np.linspace(0.0, 0.045, 70)):     # wall 6 mm to +x
        dep.add([0.006, 0, z], [0.006, 0, z + 0.0006], frame=k)

    tool = Needle.luer(inner_d=90e-6)
    tip = np.array([0.0, 0.0, 0.0])

    def worst(pose):
        return min(dep.clearance(s, frame=9999, lag=0, search=0.05)[0]
                   for _, s in pose.named_shapes())

    upright = worst(tool.at(tip))
    toward = worst(tool.at(tip, rotvec_from_axis([0.6, 0.0, 1.0])))
    away = worst(tool.at(tip, rotvec_from_axis([-0.6, 0.0, 1.0])))

    assert toward < upright, "leaning into the wall must reduce clearance"
    assert away > upright, "leaning away from it must increase clearance"


def test_from_poses_reads_six_columns(tmp_path):
    f = tmp_path / "poses.txt"
    f.write_text("\n".join([
        "# x y z rx ry rz",
        "0.0 0.0 0.000  0 0 0",
        "0.01, 0.0, 0.001,  0, 0, 0",       # commas are fine too
        "0.02\t0.0\t0.002\t0\t0\t0",        # so are tabs
        "garbage row that should be skipped",
    ]))
    p = Toolpath.from_poses(f)
    assert len(p) == 3
    assert not p.is_planar
    assert p.xyz[2] == pytest.approx([0.02, 0.0, 0.002])


def test_from_poses_scales_position_only():
    """`units` must not touch the rotation columns — radians are radians."""
    import io
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
        fh.write("10 0 5  0 0 1.5708\n")
        name = fh.name
    p = Toolpath.from_poses(name, units=1e-3)
    assert p.xyz[0] == pytest.approx([0.010, 0.0, 0.005])
    assert p.rotvec[0] == pytest.approx([0, 0, 1.5708])


# -------------------------------------------------------------------- viewer
def _collision_result():
    n = 120
    wall = np.column_stack([np.zeros(n), np.zeros(n), np.linspace(0, 0.05, n)])
    aside = np.column_stack([np.full(50, 0.003), np.zeros(50),
                             np.linspace(0.05, 0.02, 50)])
    path = Toolpath(np.vstack([wall, aside]), kinds=[PRINT] * n + [TRAVEL] * 50)
    return simulate(path, Needle.luer(inner_d=90e-6), lag=4)


def test_html_viewer_is_self_contained(tmp_path):
    """No CDN, no server — it has to open from a file:// URL years from now."""
    import re
    from toolwake import to_html
    out = to_html(_collision_result(), tmp_path / "v.html")
    h = out.read_text(encoding="utf-8")
    assert not re.findall(r"__[A-Z_]+__", h), "template placeholder left unfilled"
    assert not re.findall(r"(?:src|href)=[\"']?https?://", h), "external request"


def test_html_viewer_keeps_every_collision(tmp_path):
    """Subsampling may drop ordinary frames but never a collision.

    Losing one would make it unreachable from the buttons, which is the whole
    reason the viewer exists.
    """
    import json
    import re
    from toolwake import to_html
    res = _collision_result()
    out = to_html(res, tmp_path / "v.html", max_frames=20)   # aggressive
    D = json.loads(re.search(r"^const D = (\{.*\});$",
                             out.read_text(encoding="utf-8"), re.M).group(1))
    assert len(D["collisions"]) == res.n_hits
    for i in D["collisions"]:
        assert D["frames"][i]["c"] < 0


def test_html_worst_frame_matches_the_report(tmp_path):
    import json
    import re
    from toolwake import to_html
    res = _collision_result()
    out = to_html(res, tmp_path / "v.html")
    D = json.loads(re.search(r"^const D = (\{.*\});$",
                             out.read_text(encoding="utf-8"), re.M).group(1))
    assert D["frames"][D["worst_frame"]]["i"] == res.report()["min_clearance_row"]


def test_rpy_from_rotation_known_angles():
    from toolwake.geometry import rpy_from_rotation
    assert np.degrees(rpy_from_rotation(np.eye(3))) == pytest.approx([0, 0, 0], abs=1e-9)
    R = rotation_from_rotvec([0, 0, np.pi / 2])
    assert np.degrees(rpy_from_rotation(R)) == pytest.approx([0, 0, 90], abs=1e-6)
    R = rotation_from_rotvec([np.radians(30), 0, 0])
    assert np.degrees(rpy_from_rotation(R)) == pytest.approx([30, 0, 0], abs=1e-6)


def test_rpy_survives_gimbal_lock():
    """pitch = 90 deg makes roll and yaw the same rotation; must not blow up."""
    from toolwake.geometry import rpy_from_rotation
    r, p_, y = np.degrees(rpy_from_rotation(rotation_from_rotvec([0, np.pi / 2, 0])))
    assert p_ == pytest.approx(90, abs=1e-6)
    assert np.isfinite([r, y]).all()


def test_viewer_meta_reports_cloud_size_rpy_and_hits(tmp_path):
    import json
    import re
    from toolwake import to_html
    res = _collision_result()
    out = to_html(res, tmp_path / "v.html")
    D = json.loads(re.search(r"^const D = (\{.*\});$",
                             out.read_text(encoding="utf-8"), re.M).group(1))
    keys = [k for k, _ in D["meta"]]
    assert keys == ["point cloud size", "roll min / max", "pitch min / max",
                    "yaw min / max", "colliding frames"]
    assert dict(D["meta"])["colliding frames"] == f"{res.n_hits:,}"


def test_viewer_omits_the_part_name(tmp_path):
    """Collision rows show row + clearance only."""
    from toolwake import to_html
    h = to_html(_collision_result(), tmp_path / "v.html").read_text(encoding="utf-8")
    assert "jump to worst" not in h
    assert "COLLISION ·" not in h


def test_html_handles_a_clean_run(tmp_path):
    """No collisions must not crash the list or the navigation."""
    from toolwake import to_html
    xyz = np.column_stack([np.linspace(0, 0.05, 40), np.zeros(40), np.zeros(40)])
    res = simulate(Toolpath(xyz), Needle(inner_d=90e-6), lag=4)
    assert res.status == "ok"
    h = to_html(res, tmp_path / "clean.html").read_text(encoding="utf-8")
    assert "none — path is clear" in h


# ------------------------------------------------------------------- serving
def _loopback_carries_big_bodies(nbytes=80_000) -> bool:
    """Can this machine deliver a large response over loopback at all?

    Measured on one Windows box: 1 KB instant, 83 KB failing under both
    HTTP/1.0 and HTTP/1.1 with a bare stdlib handler and no toolwake code in
    the path — a local security product inspecting loopback, not a defect
    here. The HTTP-level test is skipped there rather than reporting a library
    failure; the same bytes are still asserted through `asgi_routes`, which is
    the integration path anyway.
    """
    import threading
    import urllib.request
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    body = b"x" * nbytes

    class H(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"

        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Connection", "close")
            self.end_headers()
            self.wfile.write(body)
            self.close_connection = True

        def log_message(self, *a):
            pass

    srv = ThreadingHTTPServer(("127.0.0.1", 0), H)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        with urllib.request.urlopen(
                f"http://127.0.0.1:{srv.server_address[1]}/", timeout=4) as r:
            return len(r.read()) == nbytes
    except Exception:
        return False
    finally:
        srv.shutdown()
        srv.server_close()


def test_view_route_returns_the_complete_page():
    """The whole page, byte for byte, with no socket in the way.

    This is the contract a host actually uses when it mounts the routes into
    its own app, so it is asserted directly rather than through TCP.
    """
    from toolwake.serve import ToolwakeSession as S, asgi_routes
    s = S(needle=Needle.luer(inner_d=90e-6))
    s.set_result(_collision_result())
    body, ctype = dict(asgi_routes(s))["/view"]()
    assert ctype.startswith("text/html")
    assert body == s.html(live=True, poll_path="version")
    assert body.rstrip().endswith("</script>")
    assert len(body) > 20_000


_LOOPBACK_OK = _loopback_carries_big_bodies()
_needs_loopback = pytest.mark.skipif(
    not _LOOPBACK_OK,
    reason="loopback here cannot carry a large body; see "
           "_loopback_carries_big_bodies. The same contract is asserted "
           "directly via asgi_routes.")



def _get(url):
    """GET a URL. HTTP status codes come back as values; transport failures
    raise, so a caller can tell a 404 apart from a dead socket."""
    import urllib.error
    import urllib.request
    try:
        with urllib.request.urlopen(url, timeout=20) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


@_needs_loopback
def test_session_serves_a_placeholder_before_any_run():
    from toolwake import ToolwakeSession
    s = ToolwakeSession()
    try:
        url = s.serve(port=0)
        code, body = _get(url + "/")
        assert code == 200 and b"No toolpath checked yet" in body
        assert _get(url + "/report.json")[0] == 404
    finally:
        s.stop()



@_needs_loopback
def test_session_serves_the_full_page_over_http():
    """End-to-end over a real socket, when the machine can manage one.

    The probe above is not sufficient on its own: the stall is intermittent, so
    a probe can succeed and the very next request still fail. A transport
    failure here says nothing about this library, so it is retried and then
    skipped — only a SHORT or MALFORMED body is treated as a real defect.
    """
    import urllib.error
    from toolwake import ToolwakeSession

    last = None
    for _ in range(3):
        s = ToolwakeSession(needle=Needle.luer(inner_d=90e-6), lag=4)
        try:
            url = s.serve(port=0)
            s.set_result(_collision_result())
            expected = len(s.html(live=True).encode("utf-8"))
            try:
                code, body = _get(url + "/")
            except (urllib.error.URLError, ConnectionError, OSError) as e:
                last = e
                continue
            assert code == 200
            assert len(body) == expected, f"served {len(body)} of {expected} bytes"
            assert body.rstrip().endswith(b"</script>")
            return
        finally:
            s.stop()
    pytest.skip(f"loopback could not carry the page after 3 tries: {last!r}")


@_needs_loopback
def test_session_version_bumps_so_an_iframe_reloads():
    import json
    from toolwake import ToolwakeSession
    s = ToolwakeSession()
    try:
        url = s.serve(port=0)
        v0 = json.loads(_get(url + "/version")[1])["version"]
        s.set_result(_collision_result())
        v1 = json.loads(_get(url + "/version")[1])["version"]
        s.set_result(_collision_result())
        v2 = json.loads(_get(url + "/version")[1])["version"]
        assert v0 < v1 < v2
    finally:
        s.stop()


def test_served_page_embeds_the_poller_but_a_file_does_not(tmp_path):
    """A static export must not poll a server that is not there."""
    from toolwake import ToolwakeSession, to_html
    s = ToolwakeSession()
    s.set_result(_collision_result())
    assert "location.reload" in s.html(live=True)
    assert "location.reload" not in s.html(live=False)
    assert "location.reload" not in to_html(
        s.result, tmp_path / "f.html").read_text(encoding="utf-8")


@_needs_loopback
def test_session_report_matches_the_result():
    import json
    from toolwake import ToolwakeSession
    s = ToolwakeSession()
    try:
        url = s.serve(port=0)
        res = _collision_result()
        s.set_result(res)
        got = json.loads(_get(url + "/report.json")[1])
        assert got["rows_penetrating"] == res.n_hits
        assert got["status"] == res.status
    finally:
        s.stop()


@_needs_loopback
def test_session_stop_is_idempotent():
    from toolwake import ToolwakeSession
    s = ToolwakeSession()
    s.serve(port=0)
    s.stop()
    s.stop()                     # must not raise
    assert s.url is None


def test_asgi_routes_shape():
    from toolwake.serve import ToolwakeSession as S, asgi_routes
    s = S()
    s.set_result(_collision_result())
    routes = dict(asgi_routes(s))
    assert set(routes) == {"/view", "/version", "/report.json"}
    body, ctype = routes["/view"]()
    assert ctype.startswith("text/html") and "location.reload" in body


# ── Broad-phase cost ────────────────────────────────────────────────────────
# These exist because a real quadratic blow-up hid behind a spatial hash that
# was working perfectly, and NOT ONE existing test noticed. Every check here
# is about how much work a query does, which is invisible to any assertion
# about what a query returns.


def test_query_examines_no_more_cells_than_necessary():
    """A query must cost min(cells in its box, cells that hold material).

    The tool's bounding box includes the luer housing, 130 mm across, while
    the grid is sized for a 45 um bead. Enumerating that box walked ~5 000
    cells per query to find a few hundred occupied ones — 43.8 million dict
    lookups over one 1 747-row toolpath, and 86% of total runtime.
    """
    d = Deposit(bead_radius=45e-6, cell=5e-3)
    for i in range(20):
        a = np.array([1e-3 * i, 0.0, 0.0])
        d.add(a, a + np.array([5e-4, 0, 0]), frame=i)

    occupied = len(d._grid)

    # A housing-sized query: vastly more cells in the box than hold material.
    big = Box(np.zeros(3), (0.065, 0.0375, 0.036))   # the 130x75x72 mm housing
    lo, hi = big.bounds(pad=0.02)
    d._candidates(lo, hi)
    box_cells = 1
    for x in (np.floor(hi / d.cell) - np.floor(lo / d.cell) + 1):
        box_cells *= int(x)
    assert box_cells > occupied, "fixture no longer exercises the large-box case"
    assert d.last_probe == occupied, (
        f"walked {d.last_probe} cells when only {occupied} hold material")

    # A needle-sized query: the box is smaller than the occupied set, so
    # enumerating it is the cheaper way round and must still be chosen.
    small = Capsule(np.array([0.0, 0.0, 0.0]), np.array([0.0, 0.0, 2e-3]), 1e-4)
    lo, hi = small.bounds(pad=1e-3)
    d._candidates(lo, hi)
    assert d.last_probe < occupied


def test_query_cost_does_not_scale_with_deposit_size():
    """Cost must follow the CANDIDATES, not the size of the whole wake.

    Storing segments in Python lists and indexing them with
    `np.asarray(self._a)[cand]` rebuilt an array of the entire deposit on
    every query. With the candidate count pinned at 2, a query went from
    673 us at 500 segments to 12.5 ms at 16 000 — pure deposit-size cost,
    none of it real work, and quadratic over a sweep.
    """
    def probe_time(n):
        d = Deposit(bead_radius=45e-6, cell=5e-3)
        # 2 cm apart, so each segment owns its own cells and a probe at the
        # origin can only ever see the first couple.
        for i in range(n):
            a = np.array([0.02 * i, 0.0, 0.0])
            d.add(a, a + np.array([1e-3, 0, 0]), frame=i)
        shape = Capsule(np.array([0.0, 0.0, 2e-3]),
                        np.array([0.0, 0.0, 1.2e-2]), 2e-4)
        d.clearance(shape, frame=n + 1, lag=10, search=0.02)   # warm up
        t0 = time.perf_counter()
        for _ in range(50):
            d.clearance(shape, frame=n + 1, lag=10, search=0.02)
        return (time.perf_counter() - t0) / 50

    small = probe_time(1000)
    large = probe_time(16000)
    # 16x the deposit for the same 2 candidates. The old code was ~11x slower
    # here; 4x leaves ample room for a loaded CI box without admitting a
    # regression of that size.
    assert large < small * 4, (
        f"query cost tracks deposit size: {small * 1e6:.0f} us at 1k segments, "
        f"{large * 1e6:.0f} us at 16k")


def test_cli_reports_its_version(capsys):
    """`toolwake --version` must work WITHOUT being given a source.

    `source` is a required positional, so this only works because argparse's
    version action exits DURING parsing, before the required-argument check
    runs. Declaration order is irrelevant — verified — but handling the flag
    after `parse_args()` instead breaks it, and the failure is confusing:

        toolwake: error: the following arguments are required: source

    which reads as a usage mistake rather than a missing feature.

    Compared against toolwake.__version__ rather than a literal, so a release
    cannot leave this test asserting last version's number.
    """
    from toolwake import __version__
    from toolwake.cli import main

    for flag in ("-v", "-V", "--version"):
        with pytest.raises(SystemExit) as e:
            main([flag])
        assert e.value.code == 0, f"{flag} should exit 0"
        assert capsys.readouterr().out.strip() == f"toolwake {__version__}"


# ── Tool solids are flat-ended ──────────────────────────────────────────────


def test_tool_tip_does_not_reach_below_itself():
    """The tip-most solid must not occupy space beneath the nozzle.

    Tool sections were `Capsule`s, and a capsule's hemispherical cap bulges a
    full radius past its endpoint. A 0.095 mm-radius cannula therefore reached
    0.095 mm BELOW its own tip — into the layer it had just printed.

    On a 0.06 mm-layer part that swallowed the layers directly underneath on
    nearly every row: 2 193 of 2 674 rows reported as collisions, worst
    -0.14 mm, every one blamed on the cannula. Material under the nozzle is
    what a print IS, so a tool model must never count it as an obstacle.
    """
    prof = Needle(inner_d=90e-6).profile
    r = prof.sections[0].r0
    layer = 60e-6
    assert layer < r, "fixture must have a layer thinner than the cannula radius"

    (_, solid), *_ = prof.capsules([0.0, 0.0, layer], [0.0, 0.0, 1.0])
    # A point one layer directly below the tip: the previous layer's bead.
    d = float(solid.distance(np.array([[0.0, 0.0, 0.0]]))[0])
    assert d > 0, (
        f"tool reports the layer {layer * 1e3:.3f} mm beneath its tip as "
        f"{d * 1e3:.4f} mm inside itself")
    # At least the vertical gap. Not exactly it: the cannula is a TUBE, so a
    # point on the axis is inside the bore and the nearest wall is the bore's
    # bottom rim — sqrt(r_inner^2 + layer^2). Asserting the gap exactly would
    # be asserting that the needle is solid.
    assert d >= layer, f"{d * 1e3:.4f} mm is less than the {layer * 1e3:.3f} mm gap"
    r_in = prof.sections[0].r_inner
    assert d == pytest.approx(np.hypot(r_in, layer), rel=1e-9)


def test_printing_above_a_previous_layer_is_not_a_collision():
    """End to end: a plain two-layer raster must not come back red."""
    layer = 60e-6
    xy = [(x * 1e-3, 0.0) for x in range(12)]
    pts, kinds = [], []
    for z in (layer, 2 * layer):
        for x, y in (xy if z == layer else xy[::-1]):
            pts.append((x, y, z))
            kinds.append(PRINT)
    path = Toolpath.from_arrays(np.array(pts), kinds=np.array(kinds))

    res = simulate(path, Needle(inner_d=90e-6), lag=2, threshold=1e-5)
    assert res.n_hits == 0, (
        f"{res.n_hits} rows called collisions on a clean two-layer raster; "
        f"worst {res.worst * 1e3:.4f} mm")


def test_bead_cannot_be_taller_than_its_layer():
    """A round bead of the bore diameter pokes up through the next layer.

    The bead defaults to the needle's bore radius, which assumes it keeps the
    cross-section it had inside the needle. Laid at a layer height smaller than
    the bore, a round bead is taller than its own layer — so its top sits above
    where the nozzle will be on the next pass, and rows that printed perfectly
    report as collisions. On a 0.036 mm-layer slab with a 0.09 mm bore that was
    498 rows, 391 of them blamed on material exactly one layer below.
    """
    layer, bore_r = 36e-6, 45e-6
    assert layer / 2 < bore_r, "fixture must have a layer thinner than the bore"

    # The second layer is offset sideways by a raster pitch, so the previous
    # layer's bead lands under the needle's WALL rather than under its bore.
    # Directly under the bore it is inside the hole and correctly ignored —
    # which is why this fixture has to be offset to reproduce anything.
    pitch = 85e-6
    assert bore_r < pitch < 95e-6, "pitch must put the bead under the wall"
    pts, kinds = [], []
    for k, z in enumerate((layer, 2 * layer)):
        for x in range(14):
            pts.append((x * 1e-3, k * pitch, z))
            kinds.append(PRINT)
    path = Toolpath.from_arrays(np.array(pts), kinds=np.array(kinds))
    needle = Needle(inner_d=2 * bore_r)

    round_bead = simulate(path, needle, lag=2, threshold=0.0)
    squashed = simulate(path, needle, lag=2, threshold=0.0,
                        bead_radius=layer / 2)

    assert round_bead.n_hits > 0, (
        "fixture no longer reproduces the artefact it exists to document")
    assert squashed.n_hits == 0, (
        f"a bead squashed to its layer still collides: {squashed.n_hits} rows, "
        f"worst {squashed.worst * 1e3:.4f} mm")
    # The clearance left is the half-layer of headroom under the nozzle.
    assert squashed.worst == pytest.approx(layer / 2, rel=1e-6)


def test_bead_sits_below_the_nozzle_face_not_on_the_path():
    """A bead centred on the toolpath is half-embedded in the nozzle's plane.

    Material leaving the bore fills the gap between the previous layer's top
    and the nozzle face: it occupies [z - h, z], centred h/2 down. Centred on
    the path instead, every bead straddles the plane the face travels in, so a
    same-layer neighbour passing under the wall reports a collision of exactly
    one bead radius — a tangent contact dressed up as a penetration. On the
    reference slab that was 41 rows, every one at exactly -bead_radius.
    """
    layer = 36e-6
    pitch = 85e-6           # neighbour lands under the wall, not the bore
    pts, kinds = [], []
    for lane in range(2):   # two adjacent lines in the SAME layer
        for x in range(14):
            pts.append((x * 1e-3, lane * pitch, layer))
            kinds.append(PRINT)
    path = Toolpath.from_arrays(np.array(pts), kinds=np.array(kinds))
    needle = Needle(inner_d=90e-6)

    on_path = simulate(path, needle, lag=2, threshold=0.0,
                       bead_radius=layer / 2)
    dropped = simulate(path, needle, lag=2, threshold=0.0,
                       bead_radius=layer / 2, bead_drop=layer / 2)

    assert on_path.n_hits > 0, "fixture no longer reproduces the artefact"
    assert on_path.worst == pytest.approx(-layer / 2, rel=1e-6), (
        "the artefact should be exactly one bead radius — a tangent contact")
    assert dropped.n_hits == 0, (
        f"bead dropped to its layer still collides: {dropped.n_hits} rows, "
        f"worst {dropped.worst * 1e3:.4f} mm")
