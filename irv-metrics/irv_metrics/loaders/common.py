"""Database helpers shared by the loaders: connection, schema, catalogue seed, load bookkeeping and cube writes.

The cube's structure comes from the YAML dataset config (contract C2) only:
`seed_catalogue(cur, catalogue)` writes sections, units, dimensions and
members, boundary sets and their levels, datasets (under their C2 ids), stored
metrics and dataset dimensions with their members, and updates rows that exist.
Every loader seeds before it writes values; `check_series` rejects a series the
catalogue does not allow.

A loader turns its source into a country-shaped object (see `CountryFacts`) and
hands it to `load_areas` and `load_facts`, inside a load record per (dataset,
boundary set, country) (`begin_load`, `finish_load`). Areas of other sets are
written with `insert_area`, related to countries with share 1 (admin, group) or
with `relate_overlay` (overlay).

DATABASE_URL selects the database (default postgresql://gri:gri@localhost:5433/gri).
"""
from __future__ import annotations

from collections import Counter
import json
import os
from pathlib import Path
import sys
from typing import Protocol

import psycopg
from psycopg.types.json import Jsonb

from ..config import Catalogue
from ..queries import CGAZ
from ..specs import CURVES

PACKAGE = Path(__file__).resolve().parent.parent
SCHEMA = PACKAGE / "schema"
DEFAULT_DSN = "postgresql://gri:gri@localhost:5433/gri"
SCHEMA_NAME = "metrics"
SEARCH_PATH = f"SET search_path = {SCHEMA_NAME}, public"


class CountryFacts(Protocol):
    """What the cube writers read from a loader's per-country result."""

    iso3: str
    boundaries: dict[str, list]                   # level -> GeoJSON features with properties adm_id, name
    names: dict[str, str]                         # adm_id -> name
    facts: dict[tuple, list]                      # (dataset, metric, dims, adm_id) -> [value, unit]
    curves: dict[tuple, list]                     # dims -> concentration curve points (specs.CURVES)


def log(*args):
    print(*args, file=sys.stderr, flush=True)


def dsn() -> str:
    return os.environ.get("DATABASE_URL", DEFAULT_DSN)


def connect() -> psycopg.Connection:
    return psycopg.connect(dsn())


ALLOW_RESET_ENV = "IRV_METRICS_ALLOW_RESET"


class SchemaManagedError(RuntimeError):
    """The host's migrations own the metrics schema; the loaders do not create or drop it."""


def schema_managed(cur) -> bool:
    """True when the database is versioned by Alembic (a `public.alembic_version` table), as in GRI.

    There the `metrics` schema is created and dropped by the host's migration
    (GRI: revision `0002_metrics`, which runs the same `schema/*.sql`), not by a loader.
    """
    cur.execute("SELECT to_regclass('public.alembic_version')")
    return cur.fetchone()[0] is not None


def apply_schema(cur, reset: bool = False):
    """Create the metrics schema if it does not exist (dropping it first with reset), then select it.

    In a database whose schema is managed by Alembic (`schema_managed`), a
    missing schema is an error (run the migrations), and `reset` is refused
    unless IRV_METRICS_ALLOW_RESET=1 (development only: drops the schema and
    recreates it from the same SQL files as the migration).
    """
    managed = schema_managed(cur)
    if reset:
        if managed and os.environ.get(ALLOW_RESET_ENV) != "1":
            raise SchemaManagedError(
                f"the {SCHEMA_NAME} schema is managed by Alembic migrations; --reset would drop it. Reload with"
                f" --if-loaded replace, or rebuild it with the migrations (alembic downgrade/upgrade), or set"
                f" {ALLOW_RESET_ENV}=1 to drop and recreate it anyway (development databases only)")
        cur.execute(f"DROP SCHEMA IF EXISTS {SCHEMA_NAME} CASCADE")
    cur.execute(f"SELECT to_regclass('{SCHEMA_NAME}.observation')")
    if cur.fetchone()[0] is None:
        if managed and not reset:
            raise SchemaManagedError(f"the {SCHEMA_NAME} schema is missing and is managed by Alembic migrations:"
                                     f" run `alembic upgrade head` first")
        log("creating schema")
        for path in sorted(SCHEMA.glob("*.sql")):
            cur.execute(path.read_text())
    cur.execute(SEARCH_PATH)


# ---------------------------------------------------------------------------
# Boundary sets and areas
# ---------------------------------------------------------------------------

def register_boundary_set(cur, boundary_set: str, spec: dict):
    """Insert or update a boundary set from a vocabulary entry (contract C2 `boundary_sets.<id>`).

    Reads `title`, `kind` (default admin), `attribution`, `attribution_url`,
    `countries.from` and `countries.min_share` (overlay; default 0.01), and
    `levels[].level`/`label` (level defaults to the entry's position).
    """
    kind = spec.get("kind", "admin")
    countries = spec.get("countries") or {}
    cur.execute("INSERT INTO boundary_set (boundary_set, kind, title, attribution, attribution_url, countries_from, min_share)"
                " VALUES (%s, %s, %s, %s, %s, %s, %s) ON CONFLICT (boundary_set) DO UPDATE SET"
                " kind = excluded.kind, title = excluded.title, attribution = excluded.attribution,"
                " attribution_url = excluded.attribution_url, countries_from = excluded.countries_from,"
                " min_share = excluded.min_share",
                (boundary_set, kind, spec["title"], spec.get("attribution"), spec.get("attribution_url"),
                 countries.get("from"), countries.get("min_share", 0.01) if kind == "overlay" else None))
    for position, level in enumerate(spec.get("levels") or []):
        cur.execute("INSERT INTO boundary_set_level VALUES (%s, %s, %s)"
                    " ON CONFLICT (boundary_set, level) DO UPDATE SET label = excluded.label",
                    (boundary_set, level.get("level", position), level.get("label")))


def insert_area(cur, boundary_set: str, area_id: str, level: int, name: str,
                geometry: dict | str | bytes | None = None, parent_area_key: int | None = None,
                source: str | None = None) -> int:
    """Insert one area and return its area_key. The geometry is GeoJSON (a dict or text) or WKB (bytes, as
    read from a GeoPackage); it is made a MultiPolygon in EPSG:4326 (WKB: polygons of a collection kept)."""
    if isinstance(geometry, (bytes, bytearray, memoryview)):
        parse = "ST_CollectionExtract(ST_GeomFromWKB(%s), 3)"
    else:
        parse = "ST_GeomFromGeoJSON(%s)"
        if geometry is not None and not isinstance(geometry, str):
            geometry = json.dumps(geometry)
    cur.execute("INSERT INTO area (boundary_set, area_id, level, name, parent_area_key, geometry, source)"
                f" VALUES (%s, %s, %s, %s, %s, ST_Multi(ST_SetSRID({parse}, 4326)), %s) RETURNING area_key",
                (boundary_set, area_id, level, name, parent_area_key, geometry, source))
    return cur.fetchone()[0]


def set_parents_by_location(cur, boundary_set: str, level: int, iso3: str | None = None):
    """Parent of each area at `level` (in one country, if given) = the nearest area at the level above
    to a point on its surface, within the same set and country."""
    cur.execute("""UPDATE area c SET parent_area_key = (
                       SELECT p.area_key FROM area p JOIN area_country pc USING (area_key)
                       WHERE p.boundary_set = c.boundary_set AND p.level = c.level - 1
                         AND pc.iso3 IN (SELECT iso3 FROM area_country WHERE area_key = c.area_key)
                       ORDER BY ST_Distance(p.geometry, ST_PointOnSurface(ST_MakeValid(c.geometry))) LIMIT 1)
                   WHERE c.boundary_set = %s AND c.level = %s AND c.geometry IS NOT NULL
                     AND (%s::text IS NULL OR c.area_key IN (SELECT area_key FROM area_country WHERE iso3 = %s))""",
                (boundary_set, level, iso3, iso3))


def relate_overlay(cur, boundary_set: str, iso3: str | None = None) -> int:
    """Relate the areas of an overlay set to the countries they overlap; return the rows written.

    Intersects each area with level 0 of the set's `countries_from` admin set and
    stores share = area(intersection) / area(area), on the spheroid, for overlaps
    of at least the set's `min_share`. With iso3, only that country is
    (re)computed, e.g. after its admin areas were loaded or reloaded; without,
    the whole set. The country rows must exist.
    """
    cur.execute("SELECT kind, countries_from, min_share FROM boundary_set WHERE boundary_set = %s", (boundary_set,))
    found = cur.fetchone()
    if found is None or found[0] != "overlay":
        raise ValueError(f"{boundary_set} is not an overlay set")
    _, admin_set, min_share = found
    cur.execute("DELETE FROM area_country ac USING area a WHERE a.area_key = ac.area_key AND a.boundary_set = %s"
                " AND (%s::text IS NULL OR ac.iso3 = %s)", (boundary_set, iso3, iso3))
    cur.execute("""
        INSERT INTO area_country (area_key, iso3, share)
        SELECT area_key, iso3, least(1.0, share) FROM (
            SELECT o.area_key, cc.iso3,
                   ST_Area(ST_Intersection(ST_MakeValid(o.geometry), ST_MakeValid(c.geometry))::geography)
                   / nullif(ST_Area(o.geometry::geography), 0) AS share
            FROM area o
            JOIN area c ON c.boundary_set = %s AND c.level = 0 AND ST_Intersects(o.geometry, c.geometry)
            JOIN area_country cc ON cc.area_key = c.area_key
            WHERE o.boundary_set = %s AND o.geometry IS NOT NULL AND (%s::text IS NULL OR cc.iso3 = %s)
        ) overlap
        WHERE share >= %s AND share > 0""", (admin_set, boundary_set, iso3, iso3, min_share))
    return cur.rowcount


def area_keys(cur, boundary_set: str, iso3: str | None = None) -> dict[str, int]:
    """area_id -> area_key for a set's areas (related to iso3, if given)."""
    cur.execute("SELECT a.area_id, a.area_key FROM area a WHERE a.boundary_set = %s"
                " AND (%s::text IS NULL OR a.area_key IN (SELECT area_key FROM area_country WHERE iso3 = %s))",
                (boundary_set, iso3, iso3))
    return dict(cur.fetchall())


def country_loaded(cur, iso3: str, boundary_set: str = CGAZ) -> bool:
    """Whether the country has areas in the (admin) set."""
    cur.execute("SELECT EXISTS (SELECT 1 FROM area_country ac JOIN area a USING (area_key)"
                " WHERE ac.iso3 = %s AND a.boundary_set = %s)", (iso3, boundary_set))
    return cur.fetchone()[0]


def remove_observations(cur, boundary_set: str, iso3: str | None = None, dataset_id: str | None = None):
    """Delete the values (and curve points) on a set's areas, optionally only those related to one
    country and of one dataset, with their load records. The scope matches a load record's key.

    Careful with overlay areas and iso3: a basin shared by KEN and UGA is related to
    both, so removing KEN's values removes the basin's values for UGA too.
    """
    scope = ("SELECT area_key FROM area WHERE boundary_set = %(set)s" if iso3 is None else
             "SELECT a.area_key FROM area a JOIN area_country ac USING (area_key)"
             " WHERE a.boundary_set = %(set)s AND ac.iso3 = %(iso3)s")
    series = "" if dataset_id is None else " AND series_id IN (SELECT series_id FROM series WHERE dataset_id = %(dataset)s)"
    params = {"set": boundary_set, "iso3": iso3, "dataset": dataset_id}
    for table in ("curve_point", "observation"):
        cur.execute(f"DELETE FROM {table} WHERE area_key IN ({scope}){series}", params)
    cur.execute("DELETE FROM load_record WHERE boundary_set = %(set)s"
                + ("" if iso3 is None else " AND iso3 = %(iso3)s")
                + ("" if dataset_id is None else " AND dataset_id = %(dataset)s"), params)


def remove_country(cur, iso3: str, boundary_set: str = CGAZ):
    """Delete a country's areas in an admin set, everything observed on them, and their load records.

    Overlay sets related to the country through `countries_from` keep their
    area_country rows; call relate_overlay(cur, set, iso3) after reloading the
    country's level 0.
    """
    cur.execute("SELECT kind FROM boundary_set WHERE boundary_set = %s", (boundary_set,))
    if (cur.fetchone() or [None])[0] != "admin":
        raise ValueError(f"{boundary_set} is not an admin set")
    remove_observations(cur, boundary_set, iso3)
    cur.execute("DELETE FROM area WHERE area_key IN (SELECT a.area_key FROM area a JOIN area_country ac USING (area_key)"
                " WHERE a.boundary_set = %s AND ac.iso3 = %s)", (boundary_set, iso3))


# ---------------------------------------------------------------------------
# Load records
# ---------------------------------------------------------------------------

def begin_load(cur, dataset_id: str, boundary_set: str, iso3: str | None, source: str, notes: dict | None = None) -> int:
    """Start the load record for (dataset, boundary set, country), replacing an earlier one; iso3 None
    for a set loaded globally. The caller removes the old values first (remove_observations)."""
    cur.execute("DELETE FROM load_record WHERE dataset_id = %s AND boundary_set = %s AND iso3 IS NOT DISTINCT FROM %s",
                (dataset_id, boundary_set, iso3))
    cur.execute("INSERT INTO load_record (dataset_id, boundary_set, iso3, source, notes) VALUES (%s, %s, %s, %s, %s)"
                " RETURNING load_id", (dataset_id, boundary_set, iso3, source, Jsonb(notes) if notes is not None else None))
    return cur.fetchone()[0]


def finish_load(cur, load_id: int):
    cur.execute("UPDATE load_record SET finished_at = clock_timestamp() WHERE load_id = %s", (load_id,))


def boundary_set_spec(catalogue: Catalogue, set_id: str) -> dict:
    """A catalogue boundary set as `register_boundary_set` reads it."""
    spec = catalogue.boundary_sets[set_id]
    return {"title": spec.title, "kind": spec.kind, "attribution": spec.attribution,
            "attribution_url": spec.attribution_url, "countries": spec.countries, "levels": list(spec.levels)}


def _upsert(cur, table: str, columns: list[str], key: list[str], rows: list[tuple]):
    """INSERT rows, updating the non-key columns of rows that exist."""
    update = [c for c in columns if c not in key]
    action = "DO UPDATE SET " + ", ".join(f"{c} = excluded.{c}" for c in update) if update else "DO NOTHING"
    cur.executemany(f"INSERT INTO {table} ({', '.join(columns)}) VALUES ({', '.join(['%s'] * len(columns))})"
                    f" ON CONFLICT ({', '.join(key)}) {action}", rows)


def seed_catalogue(cur, catalogue: Catalogue):
    """Write the catalogue's structure to the cube; rows that exist are updated to match the config.

    Seeded: sections, units, dimensions, their members (display order, labels,
    `value` as numeric_value, `code`), boundary sets and levels, datasets (C2
    ids, display order), their stored metrics (derived `share_of` metrics are
    computed by the API, not stored), dataset dimensions (default) and the
    dataset's members (with label overrides). Nothing is deleted: structure the
    config no longer has stays until the schema is rebuilt (`--reset`); it is
    logged.
    """
    _upsert(cur, "section", ["section_id", "title", "description", "sort_order"], ["section_id"],
            [(s.id, s.title, s.description, i) for i, s in enumerate(catalogue.sections.values())])
    _upsert(cur, "unit", ["unit_id", "label", "quantity"], ["unit_id"],
            [(u.id, u.label, u.quantity) for u in catalogue.units.values()])
    _upsert(cur, "dimension", ["dimension_id", "label", "description", "is_ordered"], ["dimension_id"],
            [(d.id, d.label, d.description, d.ordered) for d in catalogue.dimensions.values()])
    _upsert(cur, "dimension_member", ["dimension_id", "member_id", "label", "short_label", "description", "sort_order",
                                      "numeric_value", "code"], ["dimension_id", "member_id"],
            [(d.id, m.id, m.label, m.short_label, m.description, i, m.value, m.code)
             for d in catalogue.dimensions.values() for i, m in enumerate(d.members.values())])
    for set_id in sorted(catalogue.boundary_sets, key=lambda k: catalogue.boundary_sets[k].kind != "admin"):
        register_boundary_set(cur, set_id, boundary_set_spec(catalogue, set_id))   # admin sets first: countries_from
    datasets = list(catalogue.datasets.values())
    _upsert(cur, "dataset", ["dataset_id", "section_id", "title", "description", "documentation_path", "source_note",
                             "sort_order"], ["dataset_id"],
            [(d.id, d.section, d.title, d.description, d.documentation, d.source_note, i)
             for i, d in enumerate(datasets)])
    _upsert(cur, "metric", ["dataset_id", "metric_id", "label", "unit_id", "aggregation", "weight_metric", "description"],
            ["dataset_id", "metric_id"],
            [(d.id, m.id, m.label, m.unit, m.aggregation, m.weight, m.description)
             for d in datasets for m in d.stored_metrics.values()])
    _upsert(cur, "dataset_dimension", ["dataset_id", "dimension_id", "default_member"], ["dataset_id", "dimension_id"],
            [(d.id, k, v.default) for d in datasets for k, v in d.dimensions.items()])
    _upsert(cur, "dataset_dimension_member", ["dataset_id", "dimension_id", "member_id", "label_override"],
            ["dataset_id", "dimension_id", "member_id"],
            [(d.id, k, m, v.labels.get(m)) for d in datasets for k, v in d.dimensions.items() for m in v.members])
    cur.execute("SELECT dataset_id FROM dataset WHERE NOT dataset_id = ANY(%s) ORDER BY dataset_id",
                ([d.id for d in datasets],))
    stale = [row[0] for row in cur.fetchall()]
    if stale:
        log(f"note: datasets in the cube but not in the config (kept; --reset removes them): {stale}")


def check_series(catalogue: Catalogue, keys) -> list[str]:
    """Problems with (dataset_id, metric_id, dims as sorted pairs) series keys: the dataset, the stored
    metric, the dimensions the metric is keyed by and the dataset's members must all be in the catalogue."""
    problems = []
    for dataset_id, metric, dims in sorted(set(keys)):
        dataset = catalogue.datasets.get(dataset_id)
        if dataset is None:
            problems.append(f"{dataset_id}: not in the config")
        elif metric not in dataset.stored_metrics:
            problems.append(f"{dataset_id}: {metric} is not a stored metric")
        else:
            for dim, member in dims:
                if dim not in dataset.metrics[metric].dimensions:
                    problems.append(f"{dataset_id}: {metric} is not keyed by {dim} (={member})")
                elif member not in dataset.dimensions[dim].members:
                    problems.append(f"{dataset_id}: {dim}={member} is not one of the dataset's members")
    return problems


def load_areas(cur, country: CountryFacts, name: str, boundary_set: str = CGAZ) -> dict[str, int]:
    """Insert the country (named `name`), its admin-0 area and its boundary features into an admin
    set; relate each to the country (share 1); return area_id -> area_key for all of them."""
    iso3 = country.iso3
    cur.execute("INSERT INTO country VALUES (%s, %s, %s) ON CONFLICT (iso3) DO UPDATE SET name = excluded.name,"
                " boundary_set = coalesce(country.boundary_set, excluded.boundary_set)", (iso3, name, boundary_set))
    keys = {iso3: insert_area(cur, boundary_set, iso3, 0, name, source="union of admin-1")}
    for level in ("1", "2"):
        for feature in country.boundaries.get(level, []):
            props = feature["properties"]
            keys[props["adm_id"]] = insert_area(cur, boundary_set, props["adm_id"], int(level), props["name"],
                                                feature["geometry"], keys[iso3],  # level-2 parents set below
                                                f"GRI API /boundaries/{iso3.lower()}/{level}")
    cur.executemany("INSERT INTO area_country VALUES (%s, %s, 1)", [(key, iso3) for key in keys.values()])
    cur.execute("UPDATE area SET geometry = (SELECT ST_Multi(ST_CollectionExtract(ST_Union(ST_MakeValid(geometry)), 3))"
                " FROM area WHERE parent_area_key = %s AND level = 1) WHERE area_key = %s", (keys[iso3], keys[iso3]))
    set_parents_by_location(cur, boundary_set, 2, iso3)
    return keys


def ensure_series(cur, keys) -> dict[tuple, int]:
    """(dataset_id, metric_id, dims as sorted pairs) -> series_id, creating the missing ones in `keys`."""
    series = {}
    cur.execute("SELECT dataset_id, metric_id, dims, series_id FROM series")
    for dataset_id, metric, dims, series_id in cur.fetchall():
        series[(dataset_id, metric, tuple(sorted(dims.items())))] = series_id
    for key in sorted(set(keys) - set(series)):
        dataset_id, metric, dims = key
        cur.execute("INSERT INTO series (dataset_id, metric_id, dims) VALUES (%s, %s, %s) RETURNING series_id",
                    (dataset_id, metric, Jsonb(dict(dims))))
        series_id = cur.fetchone()[0]
        series[key] = series_id
        cur.executemany("INSERT INTO series_member VALUES (%s, %s, %s, %s)", [(series_id, dataset_id, d, m) for d, m in dims])
    return series


def fact_series(country: CountryFacts) -> set[tuple]:
    """The series keys (dataset_id, metric_id, dims) of a country's facts and curves."""
    return ({(d, m, dims) for d, m, dims, _ in country.facts}
            | {(CURVES.id, "cumulative_exposure_share", dims) for dims in country.curves})


def load_facts(cur, country: CountryFacts, areas: dict[str, int]):
    """Write a country's facts and curves; `areas` maps the facts' area ids to area keys."""
    series = ensure_series(cur, fact_series(country))
    missing = Counter()
    with cur.copy("COPY observation (series_id, area_key, value) FROM STDIN") as copy:
        for (dataset_id, metric, dims, area_id), (value, _) in country.facts.items():
            area_key = areas.get(area_id)
            if area_key is None:
                missing[area_id] += 1
                continue
            copy.write_row((series[(dataset_id, metric, dims)], area_key, value))
    if missing:
        log(f"  warning: {sum(missing.values())} facts for {len(missing)} unknown areas skipped, e.g. {list(missing)[:3]}")
    with cur.copy("COPY curve_point (series_id, area_key, seq, x, y) FROM STDIN") as copy:
        for dims, points in country.curves.items():
            for seq, point in enumerate(points):
                copy.write_row((series[(CURVES.id, "cumulative_exposure_share", dims)], areas[country.iso3], seq,
                                point["cumulative_population_share"], point["value"]))
