"""Boundary sets in the cube: admin areas and an overlay set that crosses the KEN-UGA border.

Needs a PostGIS database (DATABASE_URL); skipped without one. Everything runs in
one transaction that is rolled back, so an existing cube is left untouched (the
schema is created inside the transaction if it is missing).

Fixture geometry, simplified to rectangles that share the border at 34°E:
  KEN  34..42°E, 5°S..5°N   level 1: KEN.1 (34..38°E), KEN.2 (38..42°E)
  UGA  29..34°E, 1°S..4°N   level 1: UGA.1
  basin 1  32..36°E, 0..2°N        crosses the border: half in each country
  basin 2  33.995..38°E, 1..2°S    in Kenya but for a 0.12% sliver in Uganda,
                                   below min_share (1%), so Kenya only
"""
import os
from types import SimpleNamespace

import pytest

psycopg = pytest.importorskip("psycopg")
pytestmark = pytest.mark.skipif(not os.environ.get("DATABASE_URL"), reason="needs DATABASE_URL (PostGIS)")

from irv_metrics import queries  # noqa: E402
from irv_metrics.loaders import common  # noqa: E402

ADMIN, BASINS = "test-admin", "test-basins"


def box(x0, y0, x1, y1):
    return {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}


def feature(area_id, name, geometry):
    return {"type": "Feature", "properties": {"adm_id": area_id, "name": name}, "geometry": geometry}


class Fixture:
    def __init__(self, iso3, boundaries):
        self.iso3, self.boundaries = iso3, boundaries


@pytest.fixture(scope="module")
def db():
    with psycopg.connect(common.dsn()) as conn:
        with conn.cursor() as cur:
            common.apply_schema(cur)
            common.register_boundary_set(cur, ADMIN, {"title": "Test admin set", "levels": [{}, {}]})
            common.register_boundary_set(cur, BASINS, {"title": "Test basins", "kind": "overlay",
                                                       "countries": {"from": ADMIN, "min_share": 0.01},
                                                       "levels": [{"level": 5, "label": "Sub-basin"}]})
            for iso3, name in (("KEN", "Kenya"), ("UGA", "Uganda")):
                # Countries already in the cube keep their row; the fixture only needs them to exist.
                cur.execute("INSERT INTO country (iso3, name) VALUES (%s, %s) ON CONFLICT DO NOTHING", (iso3, name))
            common.load_areas(cur, Fixture("KEN", {"1": [feature("KEN.1", "West", box(34, -5, 38, 5)),
                                                         feature("KEN.2", "East", box(38, -5, 42, 5))]}),
                              "Kenya", ADMIN)
            common.load_areas(cur, Fixture("UGA", {"1": [feature("UGA.1", "All", box(29, -1, 34, 4))]}), "Uganda", ADMIN)
            keys = {area_id: common.insert_area(cur, BASINS, area_id, 5, f"Basin {area_id}", geometry)
                    for area_id, geometry in (("1", box(32, 0, 36, 2)), ("2", box(33.995, -2, 38, -1)))}
            related = common.relate_overlay(cur, BASINS)
            # One overlay series with a value per basin.
            cur.execute("INSERT INTO section VALUES ('other', 'Other', NULL, 99) ON CONFLICT DO NOTHING")
            cur.execute("INSERT INTO unit VALUES ('km2', 'km²', 'area') ON CONFLICT DO NOTHING")
            cur.execute("INSERT INTO dataset (dataset_id, section_id, title, sort_order) VALUES ('test/basins', 'other', 'T', 0)")
            cur.execute("INSERT INTO metric VALUES ('test/basins', 'wet_area', 'Wet area', 'km2', 'sum')")
            series_id = common.ensure_series(cur, [("test/basins", "wet_area", ())])[("test/basins", "wet_area", ())]
            load_id = common.begin_load(cur, "test/basins", BASINS, None, "fixture")
            cur.executemany("INSERT INTO observation VALUES (%s, %s, %s)",
                            [(series_id, keys["1"], 10.0), (series_id, keys["2"], 20.0)])
            common.finish_load(cur, load_id)
        yield SimpleNamespace(conn=conn, related=related, series_id=series_id)
        conn.rollback()


def test_overlay_areas_are_related_to_every_country_they_overlap(db):
    assert db.related == 3                          # basin 1 x KEN, basin 1 x UGA, basin 2 x KEN
    ken = queries.country_areas(db.conn, "KEN", BASINS, 5)
    uga = queries.country_areas(db.conn, "UGA", BASINS, 5)
    assert [a for a, _, _ in ken] == ["1", "2"]
    assert [a for a, _, _ in uga] == ["1"]                 # basin 2's sliver is below min_share
    share = {(iso3, a): s for iso3, rows in (("KEN", ken), ("UGA", uga)) for a, _, s in rows}
    assert share[("KEN", "1")] == pytest.approx(0.5, abs=1e-3)
    assert share[("UGA", "1")] == pytest.approx(0.5, abs=1e-3)
    # Areas are on the spheroid, with geodesic edges; splitting at 34°E changes them slightly.
    assert share[("KEN", "1")] + share[("UGA", "1")] == pytest.approx(1.0, abs=1e-3)
    assert share[("KEN", "2")] == pytest.approx(4 / 4.005, abs=1e-4)


def test_overlay_values_are_the_same_under_each_country(db):
    ken = queries.series_values(db.conn, db.series_id, "KEN", 5, BASINS)
    uga = queries.series_values(db.conn, db.series_id, "UGA", 5, BASINS)
    assert ken == [("1", "Basin 1", 10.0), ("2", "Basin 2", 20.0)]
    assert uga == [("1", "Basin 1", 10.0)]
    assert queries.boundaries_geojson(db.conn, "UGA", 5, BASINS).count(b'"adm_id"') == 1


def test_admin_areas_belong_to_one_country_with_share_one(db):
    assert queries.country_areas(db.conn, "KEN", ADMIN, 0) == [("KEN", "Kenya", 1.0)]
    assert queries.country_areas(db.conn, "KEN", ADMIN, 1) == [("KEN.1", "West", 1.0), ("KEN.2", "East", 1.0)]
    assert queries.country_areas(db.conn, "UGA", ADMIN, 1) == [("UGA.1", "All", 1.0)]
    assert queries.series_values(db.conn, db.series_id, "KEN", 1, ADMIN) == []     # no values on admin areas
    # Level 0 is the union of level 1; level-1 parents are the country.
    with db.conn.cursor() as cur:
        cur.execute("SELECT c.area_id, round(c.area_km2 / 1e3), p.area_id FROM area c LEFT JOIN area p"
                    " ON p.area_key = c.parent_area_key WHERE c.boundary_set = %s ORDER BY c.area_id", (ADMIN,))
        rows = cur.fetchall()
    assert [(a, p) for a, _, p in rows] == [("KEN", None), ("KEN.1", "KEN"), ("KEN.2", "KEN"), ("UGA", None), ("UGA.1", "UGA")]


def test_reloading_a_country_recomputes_its_overlaps(db):
    with db.conn.cursor() as cur:
        cur.execute("SAVEPOINT reload")
        common.relate_overlay(cur, BASINS, "UGA")
        assert {a for a, _, _ in queries.country_areas(db.conn, "UGA", BASINS, 5)} == {"1"}
        assert {a for a, _, _ in queries.country_areas(db.conn, "KEN", BASINS, 5)} == {"1", "2"}
        common.remove_observations(cur, BASINS)
        assert queries.series_values(db.conn, db.series_id, "KEN", 5, BASINS) == []
        cur.execute("SELECT count(*) FROM load_record WHERE boundary_set = %s", (BASINS,))
        assert cur.fetchone()[0] == 0
        cur.execute("ROLLBACK TO SAVEPOINT reload")
