"""Hold a result and hand it out as HTML, for embedding in another GUI.

The intended shape is an iframe in someone else's application:

    session = ToolwakeSession(needle=Needle.luer())
    url = session.serve(port=7100)            # http://127.0.0.1:7100
    ...
    session.run(toolpath)                     # iframe reloads on its own

The page polls a version endpoint, so a host that re-checks a toolpath does not
have to tell the iframe anything — it follows. That is the difference between
embedding a picture and embedding a view.

Only the standard library is used. A host already running a web framework
should skip `serve()` and hand `session.html()` to its own route instead; see
`asgi_routes` for the three lines that takes.
"""
from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from .simulate import simulate
from .tool import Needle
from .viewer import to_html_str

__all__ = ["ToolwakeSession", "asgi_routes"]


class ToolwakeSession:
    """A re-checkable simulation whose current state can be served as HTML.

    Args:
        needle: tool model reused for every run. Defaults to a plain cannula.
        title: heading shown in the viewer.
        max_frames: frames stored in the page; collisions are never dropped.
        **sim: passed to `simulate` on every `run` — `lag`, `threshold`,
            `search`.
    """

    def __init__(self, needle: Needle | None = None, *, title: str | None = None,
                 max_frames: int = 400, **sim):
        self.needle = needle or Needle()
        self.title = title
        self.max_frames = max_frames
        self.sim_kw = sim
        self._result = None
        self._version = 0
        self._lock = threading.Lock()
        self._server = None
        self._thread = None
        self._url = None

    # ------------------------------------------------------------- results
    @property
    def result(self):
        return self._result

    @property
    def version(self) -> int:
        """Bumped on every new result. The page polls this to decide to reload."""
        return self._version

    def run(self, path, **kw):
        """Simulate `path` and make the outcome the thing being served."""
        merged = dict(self.sim_kw)
        merged.update(kw)
        res = simulate(path, self.needle, **merged)
        self.set_result(res)
        return res

    def set_result(self, result) -> None:
        """Publish a result computed elsewhere."""
        with self._lock:
            self._result = result
            self._version += 1

    # --------------------------------------------------------------- views
    def html(self, *, live: bool = False, poll_path: str = "version") -> str:
        """The viewer page.

        Args:
            live: embed the poller so the page reloads when `version` changes.
                `serve()` turns this on; a static file should leave it off.
            poll_path: URL the poller hits, relative to the page.
        """
        with self._lock:
            res = self._result
        if res is None:
            return _PLACEHOLDER.replace("__TITLE__", self.title or _DEFAULT_TITLE)
        return to_html_str(res, title=self.title, max_frames=self.max_frames,
                           live_poll=poll_path if live else None)

    def report(self) -> dict | None:
        with self._lock:
            return None if self._result is None else self._result.report()

    # --------------------------------------------------------------- serve
    def serve(self, port: int = 0, host: str = "127.0.0.1") -> str:
        """Start a background HTTP server and return its URL.

        Routes: `/` and `/view` give the page, `/version` the counter the page
        polls, `/report.json` the machine-readable summary.

        `port=0` takes any free port — read the returned URL rather than
        assuming one, which keeps several sessions from fighting over a number.
        """
        if self._server is not None:
            return self._url

        session = self

        class Handler(BaseHTTPRequestHandler):
            # HTTP/1.1 gives a length-delimited body; `Connection: close` on
            # every response stops the server waiting for a follow-up request
            # that a one-shot viewer will never send.
            #
            # NOTE: on some Windows machines large localhost responses stall
            # regardless of protocol — measured here as 1 KB delivering
            # instantly while 83 KB failed under BOTH 1.0 and 1.1, with a
            # plain stdlib handler and no toolwake code involved. That is a
            # local security product inspecting loopback traffic, not a bug in
            # this server. If `serve()` hangs, mount `asgi_routes` into the
            # host application instead: same bytes, no extra socket.
            protocol_version = "HTTP/1.1"

            def _send(self, body: bytes, ctype: str, code: int = 200):
                self.send_response(code)
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(len(body)))
                # No caching: the page is regenerated per request and a cached
                # copy would defeat the whole point of the version poll.
                self.send_header("Cache-Control", "no-store")
                self.send_header("Connection", "close")
                self.end_headers()
                self.wfile.write(body)
                self.wfile.flush()
                self.close_connection = True

            def do_GET(self):
                route = self.path.split("?", 1)[0].rstrip("/") or "/"
                if route in ("/", "/view"):
                    self._send(session.html(live=True).encode("utf-8"),
                               "text/html; charset=utf-8")
                elif route == "/version":
                    self._send(json.dumps({"version": session.version}).encode(),
                               "application/json")
                elif route == "/report.json":
                    rep = session.report()
                    self._send(json.dumps(rep or {}, indent=2).encode(),
                               "application/json", 200 if rep else 404)
                else:
                    self._send(b"not found", "text/plain", 404)

            def log_message(self, *a):
                pass                    # a viewer should not spam the host's log

        self._server = ThreadingHTTPServer((host, port), Handler)
        self._url = f"http://{host}:{self._server.server_address[1]}"
        self._thread = threading.Thread(target=self._server.serve_forever,
                                        daemon=True)
        self._thread.start()
        return self._url

    @property
    def url(self) -> str | None:
        return self._url

    def stop(self) -> None:
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
            self._server = self._thread = self._url = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.stop()

    def __repr__(self):
        state = "no result" if self._result is None else repr(self._result)
        return (f"ToolwakeSession({state}, v{self._version}"
                + (f", serving {self._url}" if self._url else "") + ")")


def asgi_routes(session: ToolwakeSession):
    """Route table for mounting into an existing ASGI app, no extra port.

    Returns `[(path, handler), ...]` where each handler takes no arguments and
    returns `(body, content_type)`. Deliberately framework-agnostic so this
    module keeps its zero dependencies; wiring it to FastAPI is:

        from fastapi import Response
        for path, fn in asgi_routes(session):
            def mk(fn=fn):
                def ep():
                    body, ctype = fn()
                    return Response(body, media_type=ctype)
                return ep
            app.get("/api/toolwake" + path)(mk())

    Serving it from the host app rather than a second port avoids a separate
    process, a second port to keep open, and any cross-origin handling.
    """
    def view():
        return session.html(live=True, poll_path="version"), "text/html; charset=utf-8"

    def version():
        return json.dumps({"version": session.version}), "application/json"

    def report():
        return json.dumps(session.report() or {}, indent=2), "application/json"

    return [("/view", view), ("/version", version), ("/report.json", report)]


_DEFAULT_TITLE = "Robot arm printing needle collision checker"

_PLACEHOLDER = """<!doctype html>
<meta charset="utf-8"><title>__TITLE__</title>
<style>
  body { margin:0; height:100vh; display:flex; align-items:center;
         justify-content:center; flex-direction:column; gap:8px;
         font:14px -apple-system,Segoe UI,Helvetica,Arial,sans-serif;
         color:#6d6e71; background:#fff; }
  b { color:#1a1a1a; font-size:15px; }
</style>
<b>__TITLE__</b>
<div>No toolpath checked yet.</div>
"""
