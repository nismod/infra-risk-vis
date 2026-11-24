"""Where Alembic manages the database (GRI), the loaders neither create nor drop the metrics schema.

Needs a PostGIS database (DATABASE_URL); skipped without one. Everything runs in
one transaction that is rolled back.
"""
import os

import pytest

psycopg = pytest.importorskip("psycopg")
pytestmark = pytest.mark.skipif(not os.environ.get("DATABASE_URL"), reason="needs DATABASE_URL (PostGIS)")

from irv_metrics.loaders import common  # noqa: E402


@pytest.fixture
def cur():
    with psycopg.connect(common.dsn()) as conn:
        conn.execute("SET LOCAL lock_timeout = '5s'")
        with conn.cursor() as cur:
            try:
                yield cur
            finally:
                conn.rollback()


def managed(cur):
    cur.execute("CREATE TABLE IF NOT EXISTS public.alembic_version (version_num varchar(32) PRIMARY KEY)")


def test_unmanaged_database_creates_the_schema(cur):
    cur.execute("DROP TABLE IF EXISTS public.alembic_version")
    cur.execute("DROP SCHEMA IF EXISTS metrics CASCADE")
    assert not common.schema_managed(cur)
    common.apply_schema(cur)
    cur.execute("SELECT to_regclass('metrics.observation') IS NOT NULL")
    assert cur.fetchone()[0]


def test_managed_database_refuses_reset(cur, monkeypatch):
    monkeypatch.delenv(common.ALLOW_RESET_ENV, raising=False)
    managed(cur)
    with pytest.raises(common.SchemaManagedError, match="--reset would drop it"):
        common.apply_schema(cur, reset=True)


def test_managed_database_does_not_create_a_missing_schema(cur, monkeypatch):
    monkeypatch.delenv(common.ALLOW_RESET_ENV, raising=False)
    managed(cur)
    cur.execute("DROP SCHEMA IF EXISTS metrics CASCADE")
    with pytest.raises(common.SchemaManagedError, match="alembic upgrade head"):
        common.apply_schema(cur)


def test_managed_database_resets_with_the_dev_flag(cur, monkeypatch):
    monkeypatch.setenv(common.ALLOW_RESET_ENV, "1")
    managed(cur)
    common.apply_schema(cur, reset=True)
    cur.execute("SELECT to_regclass('metrics.observation') IS NOT NULL")
    assert cur.fetchone()[0]
