"""
The `0002_metrics` revision: the national metrics schema from irv-metrics' SQL.

Same throwaway database as test_migrations.py (DATABASE_URL); skipped without it.
"""

import importlib.util

from alembic import command
from alembic.migration import MigrationContext
from sqlalchemy import inspect, text

from test_migrations import BACKEND_DIR, autogenerate_diff, requires_db, run
from test_migrations import connection, engine  # noqa: F401  (fixtures)

# A few of the tables irv_metrics' schema/0001_init.sql creates.
METRICS_TABLES = {"area", "boundary_set", "dataset", "series", "observation", "load_record"}


def has_metrics_schema(connection) -> bool:
    return connection.execute(
        text("SELECT count(*) FROM pg_namespace WHERE nspname = 'metrics'")
    ).scalar() == 1


def test_schema_files_found_in_installed_package():
    """The revision finds the SQL through importlib.resources on the installed irv_metrics."""
    # The revision's file name is not a module name; load it by path.
    path = next((BACKEND_DIR / "migrations" / "versions").glob("*_0002_metrics.py"))
    spec = importlib.util.spec_from_file_location("rev_0002_metrics", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    names = [f.name for f in module.schema_files()]
    assert names and names == sorted(names) and "0001_init.sql" in names


@requires_db
def test_upgrade_creates_and_downgrade_drops_metrics(connection):
    run(connection, command.downgrade, "base")
    assert not has_metrics_schema(connection)

    run(connection, command.upgrade, "head")
    assert has_metrics_schema(connection)
    assert METRICS_TABLES <= set(inspect(connection).get_table_names(schema="metrics"))
    # The revision leaves the session's search_path as it found it.
    assert connection.execute(text("SELECT current_schema()")).scalar() == "public"
    # The metrics schema is excluded from autogenerate.
    assert autogenerate_diff(connection) == []

    run(connection, command.downgrade, "0001_baseline")
    assert not has_metrics_schema(connection)
    assert MigrationContext.configure(connection).get_current_revision() == "0001_baseline"

    run(connection, command.upgrade, "head")
    assert has_metrics_schema(connection)
