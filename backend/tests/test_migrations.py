"""
Migration tests against a throwaway PostGIS database.

Set DATABASE_URL to an EMPTY, disposable database (its public schema is wiped
by these tests), e.g.:

    DATABASE_URL=postgresql+psycopg2://global_dev:password@localhost:5461/global_dev pytest

Tests that need the database are skipped when DATABASE_URL is not set.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest
from alembic import command
from alembic.autogenerate import compare_metadata
from alembic.config import Config
from alembic.migration import MigrationContext
from alembic.script import ScriptDirectory
from sqlalchemy import create_engine, inspect, make_url, pool

DATABASE_URL = os.environ.get("DATABASE_URL")
BACKEND_DIR = Path(__file__).resolve().parents[1]

requires_db = pytest.mark.skipif(
    not DATABASE_URL, reason="DATABASE_URL not set (needs a throwaway PostGIS db)"
)

GRI_TABLES = {
    "adaptation_cost_benefit",
    "damages_expected",
    "damages_npv",
    "damages_rp",
    "feature_layers",
    "features",
    "raster_tile_sources",
}


@pytest.fixture(scope="module")
def engine():
    engine = create_engine(DATABASE_URL, poolclass=pool.NullPool)
    yield engine
    engine.dispose()


@pytest.fixture
def connection(engine):
    with engine.connect() as connection:
        yield connection


def alembic_config(connection) -> Config:
    config = Config(str(BACKEND_DIR / "alembic.ini"))
    config.attributes["connection"] = connection
    config.attributes["configure_logger"] = False
    return config


def run(connection, fn, *args):
    fn(alembic_config(connection), *args)
    connection.commit()


def autogenerate_diff(connection):
    """What `alembic revision --autogenerate` would emit, using env.py's filters."""
    import importlib.util

    from alembic.runtime.environment import EnvironmentContext

    config = alembic_config(connection)
    script = ScriptDirectory.from_config(config)
    diff = []

    def compare(rev, context):
        diff.extend(compare_metadata(context, context.opts["target_metadata"]))
        return []

    with EnvironmentContext(config, script, fn=compare):
        spec = importlib.util.spec_from_file_location(
            "migrations_env", BACKEND_DIR / "migrations" / "env.py"
        )
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    return diff


def gri_tables(connection):
    return set(inspect(connection).get_table_names()) & GRI_TABLES


def test_single_head():
    script = ScriptDirectory.from_config(Config(str(BACKEND_DIR / "alembic.ini")))
    assert len(script.get_heads()) == 1


@requires_db
def test_upgrade_then_autogenerate_is_empty(connection):
    run(connection, command.downgrade, "base")
    assert gri_tables(connection) == set()

    run(connection, command.upgrade, "head")
    assert gri_tables(connection) == GRI_TABLES
    current = MigrationContext.configure(connection).get_current_revision()
    script = ScriptDirectory.from_config(alembic_config(connection))
    assert current == script.get_current_head()

    assert autogenerate_diff(connection) == []


@requires_db
def test_downgrade_base_then_upgrade_head(connection):
    run(connection, command.upgrade, "head")
    run(connection, command.downgrade, "base")
    assert gri_tables(connection) == set()
    assert MigrationContext.configure(connection).get_current_revision() is None

    run(connection, command.upgrade, "head")
    assert gri_tables(connection) == GRI_TABLES
    assert autogenerate_diff(connection) == []


def alembic_cli(*args):
    """Run the alembic CLI as the migrate service does: PG* env vars, no URL."""
    url = make_url(DATABASE_URL)
    env = {
        **os.environ,
        "PGHOST": url.host or "localhost",
        "PGPORT": str(url.port or 5432),
        "PGUSER": url.username or "",
        "PGPASSWORD": url.password or "",
        "PGDATABASE": url.database or "",
    }
    return subprocess.run(
        [sys.executable, "-m", "alembic", *args],
        cwd=BACKEND_DIR,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )


@requires_db
def test_cli_upgrade_commits(engine):
    """`alembic downgrade base` / `upgrade head` via the CLI persist their changes."""
    alembic_cli("downgrade", "base")
    with engine.connect() as connection:
        assert gri_tables(connection) == set()

    alembic_cli("upgrade", "head")
    with engine.connect() as connection:
        assert gri_tables(connection) == GRI_TABLES

    assert "No new upgrade operations detected" in alembic_cli("check").stdout
