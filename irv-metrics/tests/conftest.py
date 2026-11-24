"""Shared fixtures.

`api_base`: the URL of a running dev app (irv_metrics.adapters.fastapi). Uses
IRV_METRICS_BASE_URL when set; otherwise starts the app in-process on a free
port, against DATABASE_URL. Skips without DATABASE_URL, or when the database
has no loaded country (load one with `python -m irv_metrics.loaders.scrape`).
"""
import os
import socket
import threading
import time
from urllib.request import urlopen

import pytest


def _free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return sock.getsockname()[1]


def serve(app):
    """Run an ASGI app with uvicorn in a daemon thread; return (base URL, server)."""
    import uvicorn

    port = _free_port()
    server = uvicorn.Server(uvicorn.Config(app, host="127.0.0.1", port=port, log_level="warning"))
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    deadline = time.monotonic() + 30
    while not server.started:
        if time.monotonic() > deadline or not thread.is_alive():
            pytest.fail("the API did not start")
        time.sleep(0.05)
    server.thread = thread
    return f"http://127.0.0.1:{port}", server


def stop(server):
    server.should_exit = True
    server.thread.join(10)


@pytest.fixture(scope="session")
def serve_app():
    """serve_app(app) -> base URL; the servers stop at the end of the session."""
    servers = []

    def start(app):
        base, server = serve(app)
        servers.append(server)
        return base

    yield start
    for server in servers:
        stop(server)


@pytest.fixture(scope="session")
def api_base():
    base = os.environ.get("IRV_METRICS_BASE_URL")
    if base:
        yield base.rstrip("/")
        return
    if not os.environ.get("DATABASE_URL"):
        pytest.skip("needs DATABASE_URL (a loaded cube) or IRV_METRICS_BASE_URL")
    import psycopg

    try:
        with psycopg.connect(os.environ["DATABASE_URL"]) as conn:
            loaded = conn.execute("SELECT to_regclass('metrics.load_record') IS NOT NULL").fetchone()[0] and \
                conn.execute("SELECT count(*) FROM metrics.load_record WHERE finished_at IS NOT NULL").fetchone()[0]
    except psycopg.OperationalError as error:
        pytest.skip(f"database not reachable: {error}")
    if not loaded:
        pytest.skip("the database has no loaded data")

    from irv_metrics.adapters.fastapi import app

    base, server = serve(app)
    with urlopen(base + "/health", timeout=10):
        pass
    yield base
    stop(server)
