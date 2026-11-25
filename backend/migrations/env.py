"""
Alembic environment for the GRI backend database.

To connect to the database, set the ``PG*`` environment variables (see
https://www.postgresql.org/docs/current/libpq-envars.html).

Tests may instead pass an open connection as
``config.attributes["connection"]``.

Autogenerate compares ``backend.db.models.Base.metadata`` against the default
schema only. Objects that belong to database extensions (PostGIS's
``spatial_ref_sys``, ``pg_cron``'s tables, ...) and schemas that are managed
elsewhere are filtered out, so they never appear in a diff.
"""

from logging.config import fileConfig

from alembic import context
from geoalchemy2 import alembic_helpers
from sqlalchemy import create_engine, pool, text
from sqlalchemy.engine import Connection

from backend.db import models
from backend.db.database import engine as app_engine

config = context.config

if config.config_file_name is not None and config.attributes.get(
    "configure_logger", True
):
    fileConfig(config.config_file_name)

target_metadata = models.Base.metadata

# Schemas not managed by autogenerate. ``include_schemas`` is off, so only the
# default schema is compared anyway for now.
# - topology, tiger, tiger_data: PostGIS extensions
# - cron: pg_cron (installed in production)
# - metrics: national metrics cube, created from irv-metrics' SQL files by
#   hand-written revisions, not from models.py
EXCLUDED_SCHEMAS = {"topology", "tiger", "tiger_data", "cron", "metrics"}

# Extension tables excluded by name, in case they are not registered as
# extension members (e.g. restored from a dump without CREATE EXTENSION).
EXCLUDED_TABLES = {
    "spatial_ref_sys",
    "geometry_columns",
    "geography_columns",
    "raster_columns",
    "raster_overviews",
}

EXTENSION_OWNED_TABLES_SQL = text(
    """
    SELECT n.nspname, c.relname, pg_catalog.pg_table_is_visible(c.oid)
    FROM pg_depend d
    JOIN pg_class c ON d.classid = 'pg_class'::regclass AND d.objid = c.oid
    JOIN pg_namespace n ON n.oid = c.relnamespace
    WHERE d.deptype = 'e' AND d.refclassid = 'pg_extension'::regclass
    """
)


def _make_filters(connection: Connection):
    """Return the ``include_name`` filter for autogenerate."""
    # (schema, table) of extension member tables, plus the names of those
    # visible on the search_path.
    #
    # Reflection with schema=None lists every visible table (e.g. topology.layer
    # when the search_path includes topology), and reports them without a
    # schema, so match those by name.
    extension_tables: tuple[set[tuple[str, str]], set[str]] | None = None

    def is_extension_table(schema: str | None, name: str) -> bool:
        # Queried lazily, i.e. only during autogenerate and inside the
        # migration transaction. Querying before context.begin_transaction()
        # would autobegin a transaction that Alembic then never commits.
        nonlocal extension_tables
        if extension_tables is None:
            rows = connection.execute(EXTENSION_OWNED_TABLES_SQL).all()
            extension_tables = (
                {(row[0], row[1]) for row in rows},
                {row[1] for row in rows if row[2]},
            )
        qualified, visible = extension_tables
        if schema is None:
            return name in visible
        return (schema, name) in qualified

    def include_name(name, type_, parent_names):
        if type_ == "schema":
            return name not in EXCLUDED_SCHEMAS
        if type_ == "table":
            schema = parent_names.get("schema_name")
            return (
                schema not in EXCLUDED_SCHEMAS
                and name not in EXCLUDED_TABLES
                and not is_extension_table(schema, name)
            )
        return True

    return include_name


def configure_context(connection: Connection) -> None:
    context.configure(
        connection=connection,
        target_metadata=target_metadata,
        include_name=_make_filters(connection),
        # geoalchemy2's filter for spatial extension tables, by name
        include_object=alembic_helpers.include_object,
        render_item=alembic_helpers.render_item,
        compare_type=True,
        compare_server_default=True,
    )


def run_migrations_offline() -> None:
    """Emit SQL to stdout instead of running it (``alembic upgrade --sql``)."""
    context.configure(
        url=app_engine.url,
        target_metadata=target_metadata,
        literal_binds=True,
        dialect_opts={"paramstyle": "named"},
        render_item=alembic_helpers.render_item,
    )
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    connection = config.attributes.get("connection")
    if connection is not None:
        configure_context(connection)
        with context.begin_transaction():
            context.run_migrations()
        return

    connectable = create_engine(app_engine.url, poolclass=pool.NullPool)
    with connectable.connect() as connection:
        configure_context(connection)
        with context.begin_transaction():
            context.run_migrations()


if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()
