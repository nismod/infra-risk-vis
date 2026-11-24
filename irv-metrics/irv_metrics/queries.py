"""Every SQL query the API runs, as plain functions over a psycopg connection.

The connection must have search_path = metrics, public. No web framework here:
functions take arguments and return rows or plain Python data.

Areas are named by (boundary set, area_id); country-scoped queries filter
through area_country, so an overlay area appears under every country it
overlaps. Functions that serve the upstream-compatible routes default to the
CGAZ admin set, which is what the scraped data is keyed on.
"""
from __future__ import annotations

import json

import psycopg

CGAZ = "cgaz-6.0.0"

# The country scope: areas of one set at one level related to one country, through area_country.
# Use as "FROM area_country ac JOIN area a USING (area_key) ... WHERE {IN_COUNTRY}".
IN_COUNTRY = "ac.iso3 = %(iso3)s AND a.boundary_set = %(set)s AND a.level = %(level)s"


# ---------------------------------------------------------------------------
# Catalog: small, rarely changing tables the API keeps in memory
# ---------------------------------------------------------------------------

def load_version(conn: psycopg.Connection) -> tuple:
    """(latest finished load, number of finished loads, latest finish time): changes whenever a load completes."""
    return conn.execute("SELECT max(load_id), count(*), max(finished_at) FROM load_record"
                        " WHERE finished_at IS NOT NULL").fetchone()


def series_keys(conn: psycopg.Connection) -> dict[tuple, int]:
    """(dataset_id, metric_id, frozenset(dims.items())) -> series_id."""
    return {(d, m, frozenset(dims.items())): s
            for s, d, m, dims in conn.execute("SELECT series_id, dataset_id, metric_id, dims FROM series")}


def metric_units(conn: psycopg.Connection) -> dict[tuple, str]:
    """(dataset_id, metric_id) -> unit_id."""
    return {(d, m): u for d, m, u in conn.execute("SELECT dataset_id, metric_id, unit_id FROM metric")}


def admin_areas(conn: psycopg.Connection, boundary_set: str = CGAZ) -> dict[str, tuple[str, int, str]]:
    """area_id -> (iso3, level, name) for the areas of an admin set (each has one country)."""
    return {a: (iso3.strip(), level, name) for a, iso3, level, name in conn.execute(
        "SELECT a.area_id, ac.iso3, a.level, a.name FROM area a JOIN area_country ac USING (area_key)"
        " WHERE a.boundary_set = %s", (boundary_set,))}


def member_order(conn: psycopg.Connection) -> dict[tuple, int]:
    """(dimension_id, member_id) -> sort_order."""
    return {(d, m): o for d, m, o in conn.execute("SELECT dimension_id, member_id, sort_order FROM dimension_member")}


def countries(conn: psycopg.Connection) -> list[dict]:
    return [{"iso3": i.strip(), "name": n} for i, n in conn.execute("SELECT iso3, name FROM country ORDER BY iso3")]


def country_areas(conn: psycopg.Connection, iso3: str, boundary_set: str, level: int) -> list[tuple]:
    """(area_id, name, share) for the areas of a set at one level related to a country, by area_id."""
    return conn.execute("SELECT a.area_id, a.name, ac.share FROM area_country ac JOIN area a USING (area_key)"
                        f" WHERE {IN_COUNTRY} ORDER BY a.area_id",
                        {"iso3": iso3, "set": boundary_set, "level": level}).fetchall()


# ---------------------------------------------------------------------------
# Canonical API
# ---------------------------------------------------------------------------

def dataset_structure(conn: psycopg.Connection) -> list[dict]:
    """Every dataset with its metrics (and units) and the members of each dimension."""
    rows = conn.execute("""
        SELECT d.dataset_id, d.section_id, d.title,
               (SELECT jsonb_agg(jsonb_build_object('id', m.metric_id, 'label', m.label, 'unit', m.unit_id)
                                 ORDER BY m.metric_id) FROM metric m WHERE m.dataset_id = d.dataset_id),
               (SELECT jsonb_object_agg(dd.dimension_id, (
                    SELECT jsonb_agg(jsonb_build_object('id', dm.member_id, 'label', coalesce(ddm.label_override, dm.label))
                                     ORDER BY dm.sort_order)
                    FROM dataset_dimension_member ddm JOIN dimension_member dm USING (dimension_id, member_id)
                    WHERE ddm.dataset_id = d.dataset_id AND ddm.dimension_id = dd.dimension_id))
                FROM dataset_dimension dd WHERE dd.dataset_id = d.dataset_id)
        FROM dataset d ORDER BY d.sort_order""").fetchall()
    return [{"id": d, "section": s, "title": t, "metrics": m or [], "dimensions": dims or {}} for d, s, t, m, dims in rows]


def observations(conn: psycopg.Connection, dataset: str, iso3: str, *, level: int | None = None,
                 adm_id: str | None = None, metric: str | None = None, filters: dict | None = None,
                 limit: int = 100_000, boundary_set: str = CGAZ) -> list[dict]:
    """Long-format values of one dataset in one country; `filters` are {dimension: member}."""
    sql = ["SELECT a.area_id, a.name, a.level, s.metric_id, s.dims, o.value, m.unit_id",
           "FROM observation o JOIN series s USING (series_id) JOIN area a USING (area_key)",
           "JOIN area_country ac USING (area_key) JOIN metric m USING (dataset_id, metric_id)",
           "WHERE s.dataset_id = %s AND ac.iso3 = %s AND a.boundary_set = %s"]
    params: list = [dataset, iso3, boundary_set]
    for clause, value in (("a.level = %s", level), ("a.area_id = %s", adm_id), ("s.metric_id = %s", metric)):
        if value is not None:
            sql.append(f"AND {clause}")
            params.append(value)
    if filters:
        sql.append("AND s.dims @> %s::jsonb")
        params.append(json.dumps(filters))
    sql.append("ORDER BY s.series_id, a.area_id LIMIT %s")
    params.append(limit)
    rows = conn.execute(" ".join(sql), params).fetchall()
    return [{"adm_id": a, "name": n, "level": lv, "metric": m, "dims": d, "value": v, "unit": u}
            for a, n, lv, m, d, v, u in rows]


# ---------------------------------------------------------------------------
# Upstream-compatible routes
# ---------------------------------------------------------------------------

def series_values(conn: psycopg.Connection, series_id: int, iso3: str, level: int,
                  boundary_set: str = CGAZ) -> list[tuple]:
    """(area_id, name, value) for one series over the areas of a set at one level related to a country."""
    return conn.execute(
        "SELECT a.area_id, a.name, o.value FROM area_country ac JOIN area a USING (area_key)"
        f" JOIN observation o USING (area_key) WHERE {IN_COUNTRY} AND o.series_id = %(series)s ORDER BY a.area_id",
        {"series": series_id, "iso3": iso3, "set": boundary_set, "level": level}).fetchall()


def area_observations(conn: psycopg.Connection, dataset_id: str, area_id: str,
                      boundary_set: str = CGAZ) -> list[tuple]:
    """(metric_id, dims, value) for every observation of one dataset on one area."""
    return conn.execute("SELECT s.metric_id, s.dims, o.value FROM area a JOIN observation o USING (area_key)"
                        " JOIN series s USING (series_id)"
                        " WHERE s.dataset_id = %s AND a.boundary_set = %s AND a.area_id = %s",
                        (dataset_id, boundary_set, area_id)).fetchall()


def download_rows(conn: psycopg.Connection, dataset_id: str, iso3: str, level: int,
                  boundary_set: str = CGAZ) -> list[tuple]:
    """(area_id, name, metric_id, dims, unit_id, value) for a card CSV, by area then series."""
    return conn.execute(
        "SELECT a.area_id, a.name, s.metric_id, s.dims, m.unit_id, o.value FROM area_country ac JOIN area a USING (area_key)"
        " JOIN observation o USING (area_key) JOIN series s USING (series_id) JOIN metric m USING (dataset_id, metric_id)"
        f" WHERE {IN_COUNTRY} AND s.dataset_id = %(dataset)s ORDER BY a.area_id, s.series_id",
        {"dataset": dataset_id, "iso3": iso3, "set": boundary_set, "level": level}).fetchall()


def boundaries_geojson(conn: psycopg.Connection, iso3: str, level: int, boundary_set: str = CGAZ) -> bytes | None:
    """GeoJSON FeatureCollection (encoded) of the areas of a set at one level related to a country;
    None when there are none."""
    rows = conn.execute("SELECT a.area_id, a.name, ST_AsGeoJSON(a.geometry, 6) FROM area_country ac"
                        f" JOIN area a USING (area_key) WHERE {IN_COUNTRY} ORDER BY a.area_id", {"iso3": iso3, "set": boundary_set, "level": level}).fetchall()
    if not rows:
        return None
    features = ",".join(
        '{"type":"Feature","id":%s,"properties":{"adm_id":%s,"name":%s},"geometry":%s}'
        % (json.dumps(a), json.dumps(a), json.dumps(n), g) for a, n, g in rows)
    return ('{"type":"FeatureCollection","features":[%s]}' % features).encode()


def curve_points(conn: psycopg.Connection, series_id: int, area_id: str, boundary_set: str = CGAZ) -> list[tuple]:
    """(x, y) points of one concentration curve, in order."""
    return conn.execute("SELECT p.x, p.y FROM curve_point p JOIN area a USING (area_key)"
                        " WHERE p.series_id = %s AND a.boundary_set = %s AND a.area_id = %s ORDER BY p.seq",
                        (series_id, boundary_set, area_id)).fetchall()


# ---------------------------------------------------------------------------
# /v1 (contract C1)
#
# Areas are addressed by (boundary_set, area_id) and scoped to a country through
# area_country; nothing here assumes the CGAZ set.
# ---------------------------------------------------------------------------

def series_list(conn: psycopg.Connection) -> list[tuple]:
    """(series_id, dataset_id, metric_id, dims) for every series."""
    return conn.execute("SELECT series_id, dataset_id, metric_id, dims FROM series ORDER BY series_id").fetchall()


def boundary_set_rows(conn: psycopg.Connection) -> list[tuple]:
    """(boundary_set, kind, title) for every set."""
    return conn.execute("SELECT boundary_set, kind, title FROM boundary_set ORDER BY boundary_set").fetchall()


def boundary_set_level_labels(conn: psycopg.Connection) -> dict[tuple[str, int], str | None]:
    """(boundary_set, level) -> label."""
    return {(s, lv): label for s, lv, label in conn.execute("SELECT boundary_set, level, label FROM boundary_set_level")}


def country_rows(conn: psycopg.Connection) -> list[tuple]:
    """(iso3, name, admin set, west, south, east, north) for every country. The box is the extent of the
    country's level-0 area in its admin set, or of all its areas in that set when level 0 has no geometry."""
    return [(iso3.strip(), *rest) for iso3, *rest in conn.execute("""
        SELECT c.iso3, c.name, c.boundary_set,
               ST_XMin(e), ST_YMin(e), ST_XMax(e), ST_YMax(e)
        FROM country c
        LEFT JOIN LATERAL (
            SELECT coalesce(
                (SELECT ST_Extent(a.geometry) FROM area a JOIN area_country ac USING (area_key)
                 WHERE ac.iso3 = c.iso3 AND a.boundary_set = c.boundary_set AND a.level = 0),
                (SELECT ST_Extent(a.geometry) FROM area a JOIN area_country ac USING (area_key)
                 WHERE ac.iso3 = c.iso3 AND a.boundary_set = c.boundary_set))::geometry AS e
        ) extent ON true
        ORDER BY c.iso3""")]


def country_set_levels(conn: psycopg.Connection) -> list[tuple]:
    """(iso3, boundary_set, level, number of areas) for every set and level with areas related to a country."""
    return [(iso3.strip(), s, lv, n) for iso3, s, lv, n in conn.execute(
        "SELECT ac.iso3, a.boundary_set, a.level, count(*) FROM area_country ac JOIN area a USING (area_key)"
        " GROUP BY 1, 2, 3 ORDER BY 1, 2, 3")]


def coverage(conn: psycopg.Connection) -> list[tuple]:
    """(iso3, dataset_id, boundary_set, level) wherever a dataset has values (or curve points) on the
    areas of a set at a level related to a country. Scans the facts once; the API runs it after each load."""
    return [(iso3.strip(), d, s, lv) for iso3, d, s, lv in conn.execute("""
        WITH present AS (
            SELECT DISTINCT o.series_id, a.boundary_set, a.level, o.area_key
            FROM (SELECT series_id, area_key FROM observation UNION ALL
                  SELECT DISTINCT series_id, area_key FROM curve_point) o
            JOIN area a USING (area_key))
        SELECT DISTINCT ac.iso3, s.dataset_id, p.boundary_set, p.level
        FROM present p JOIN series s USING (series_id) JOIN area_country ac USING (area_key)""")]


def areas_in_country(conn: psycopg.Connection, boundary_set: str, area_ids: list[str], iso3: str) -> list[tuple]:
    """(area_id, area_key, level, name, share) for those of `area_ids` in the set that are related to the country."""
    return conn.execute(
        "SELECT a.area_id, a.area_key, a.level, a.name, ac.share FROM area a JOIN area_country ac USING (area_key)"
        " WHERE a.boundary_set = %s AND a.area_id = ANY(%s) AND ac.iso3 = %s",
        (boundary_set, list(area_ids), iso3)).fetchall()


def series_area_values(conn: psycopg.Connection, series_ids: list[int], iso3: str, boundary_set: str,
                       level: int) -> list[tuple]:
    """(area_id, name, share, series_id, value) for the areas of a set at one level related to a country,
    by area_id, with one row per area and series in `series_ids`. Value is None where nothing is stored."""
    return conn.execute(
        "SELECT a.area_id, a.name, ac.share, s.series_id, o.value"
        " FROM area_country ac JOIN area a USING (area_key) CROSS JOIN unnest(%(series)s::int4[]) s(series_id)"
        " LEFT JOIN observation o ON o.area_key = a.area_key AND o.series_id = s.series_id"
        f" WHERE {IN_COUNTRY} ORDER BY a.area_id",
        {"series": list(series_ids), "iso3": iso3, "set": boundary_set, "level": level}).fetchall()


def area_series_values(conn: psycopg.Connection, area_keys: list[int], series_ids: list[int]) -> list[tuple]:
    """(area_key, series_id, value) for every stored value of the given series on the given areas."""
    return conn.execute("SELECT area_key, series_id, value FROM observation"
                        " WHERE area_key = ANY(%s) AND series_id = ANY(%s)",
                        (list(area_keys), list(series_ids))).fetchall()


def curve(conn: psycopg.Connection, series_id: int, area_key: int) -> list[tuple]:
    """(x, y) points of one curve on one area, in order."""
    return conn.execute("SELECT x, y FROM curve_point WHERE series_id = %s AND area_key = %s ORDER BY seq",
                        (series_id, area_key)).fetchall()


def level_observations(conn: psycopg.Connection, series_ids: list[int], iso3: str, boundary_set: str,
                       level: int) -> list[tuple]:
    """(area_id, name, share, series_id, value) for every stored value of the given series on the areas
    of a set at one level related to a country, by area_id."""
    return conn.execute(
        "SELECT a.area_id, a.name, ac.share, o.series_id, o.value FROM area_country ac JOIN area a USING (area_key)"
        f" JOIN observation o USING (area_key) WHERE {IN_COUNTRY} AND o.series_id = ANY(%(series)s)"
        " ORDER BY a.area_id",
        {"series": list(series_ids), "iso3": iso3, "set": boundary_set, "level": level}).fetchall()


def level_curve_points(conn: psycopg.Connection, series_ids: list[int], iso3: str, boundary_set: str,
                       level: int) -> list[tuple]:
    """(area_id, name, share, series_id, x, y) for every curve point of the given series on the areas of a
    set at one level related to a country, by area_id then point order."""
    return conn.execute(
        "SELECT a.area_id, a.name, ac.share, p.series_id, p.x, p.y FROM area_country ac JOIN area a USING (area_key)"
        f" JOIN curve_point p USING (area_key) WHERE {IN_COUNTRY} AND p.series_id = ANY(%(series)s)"
        " ORDER BY a.area_id, p.series_id, p.seq",
        {"series": list(series_ids), "iso3": iso3, "set": boundary_set, "level": level}).fetchall()


def area_features(conn: psycopg.Connection, iso3: str, boundary_set: str, level: int,
                  with_share: bool) -> bytes | None:
    """GeoJSON FeatureCollection (encoded) of the areas of a set at one level related to a country, by
    area_id: `id` = area_id; properties area_id, name, parent_area_id, and share when `with_share`.
    Whole geometries, WGS84, 6 decimal places. None when there are no areas."""
    rows = conn.execute(
        "SELECT a.area_id, a.name, p.area_id, ac.share, ST_AsGeoJSON(a.geometry, 6)"
        " FROM area_country ac JOIN area a USING (area_key) LEFT JOIN area p ON p.area_key = a.parent_area_key"
        f" WHERE {IN_COUNTRY} ORDER BY a.area_id", {"iso3": iso3, "set": boundary_set, "level": level}).fetchall()
    if not rows:
        return None
    features = []
    for area_id, name, parent, share, geometry in rows:
        props = {"area_id": area_id, "name": name, "parent_area_id": parent}
        if with_share:
            props["share"] = share
        features.append('{"type":"Feature","id":%s,"properties":%s,"geometry":%s}'
                        % (json.dumps(area_id), json.dumps(props, ensure_ascii=False), geometry or "null"))
    return ('{"type":"FeatureCollection","features":[%s]}' % ",".join(features)).encode()
