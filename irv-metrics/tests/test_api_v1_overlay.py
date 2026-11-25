"""The /v1 service on an overlay set and a group set (the scrape has only CGAZ).

Fixture as in test_boundary_sets.py: two admin countries split at 34°E and two
basins, basin 1 half in each country, basin 2 in Kenya only; plus a group set
with one group of both countries. A hand-built config catalogue has one dataset
on the basins. Everything runs in one transaction that is rolled back; the
service reads through that connection.

Needs a PostGIS database (DATABASE_URL); skipped without one.
"""
from contextlib import contextmanager
import csv
import io
import json
import os

import pytest

psycopg = pytest.importorskip("psycopg")
pytestmark = pytest.mark.skipif(not os.environ.get("DATABASE_URL"), reason="needs DATABASE_URL (PostGIS)")

from irv_metrics import config  # noqa: E402
from irv_metrics.api import ApiError, Service  # noqa: E402
from irv_metrics.loaders import common  # noqa: E402

ADMIN, BASINS, GROUPS = "test-admin", "test-basins", "test-groups"
DATASET = "test/basins"


def box(x0, y0, x1, y1):
    return {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}


def feature(area_id, name, geometry):
    return {"type": "Feature", "properties": {"adm_id": area_id, "name": name}, "geometry": geometry}


class Fixture:
    def __init__(self, iso3, boundaries):
        self.iso3, self.boundaries = iso3, boundaries


class OneConnection:
    """A pool that hands out one connection (the test's transaction)."""

    def __init__(self, conn):
        self.conn = conn

    @contextmanager
    def connection(self):
        yield self.conn


def catalogue() -> config.Catalogue:
    sets = {
        ADMIN: config.BoundarySet(ADMIN, "Test admin set", "admin", ({}, {})),
        BASINS: config.BoundarySet(BASINS, "Test basins", "overlay", ({"level": 5, "label": "Sub-basin"},),
                                   {"from": ADMIN, "min_share": 0.01}),
        GROUPS: config.BoundarySet(GROUPS, "Test groups", "group", ({"level": 0, "label": "Group"},)),
    }
    dataset = config.Dataset(DATASET, "Basins", {"wet_area": config.Metric("wet_area", "Wet area", "km2", (), "sum")},
                             {}, sources=({"boundary_set": BASINS, "files": "x", "format": "long"},
                                          {"boundary_set": GROUPS, "files": "y", "format": "long"}))
    return config.Catalogue({"test": config.Section("test", "Test")}, {"km2": config.Unit("km2", "km²", "area")},
                            {}, sets, {}, {DATASET: dataset})


@pytest.fixture(scope="module")
def service():
    with psycopg.connect(common.dsn()) as conn:
        with conn.cursor() as cur:
            common.apply_schema(cur)
            common.register_boundary_set(cur, ADMIN, {"title": "Test admin set", "levels": [{}, {}]})
            common.register_boundary_set(cur, BASINS, {"title": "Test basins", "kind": "overlay",
                                                       "countries": {"from": ADMIN, "min_share": 0.01},
                                                       "levels": [{"level": 5, "label": "Sub-basin"}]})
            common.register_boundary_set(cur, GROUPS, {"title": "Test groups", "kind": "group",
                                                       "levels": [{"label": "Group"}]})
            for iso3, name in (("KEN", "Kenya"), ("UGA", "Uganda")):
                cur.execute("INSERT INTO country (iso3, name) VALUES (%s, %s) ON CONFLICT DO NOTHING", (iso3, name))
            common.load_areas(cur, Fixture("KEN", {"1": [feature("KEN.1", "West", box(34, -5, 38, 5)),
                                                         feature("KEN.2", "East", box(38, -5, 42, 5))]}),
                              "Kenya", ADMIN)
            common.load_areas(cur, Fixture("UGA", {"1": [feature("UGA.1", "All", box(29, -1, 34, 4))]}), "Uganda", ADMIN)
            keys = {area_id: common.insert_area(cur, BASINS, area_id, 5, f"Basin {area_id}", geometry)
                    for area_id, geometry in (("1", box(32, 0, 36, 2)), ("2", box(33.995, -2, 38, -1)))}
            common.relate_overlay(cur, BASINS)
            group = common.insert_area(cur, GROUPS, "EAST", 0, "East Africa")
            cur.executemany("INSERT INTO area_country VALUES (%s, %s, 1)", [(group, "KEN"), (group, "UGA")])
            cur.execute("INSERT INTO section VALUES ('other', 'Other', NULL, 99) ON CONFLICT DO NOTHING")
            cur.execute("INSERT INTO unit VALUES ('km2', 'km²', 'area') ON CONFLICT DO NOTHING")
            cur.execute("INSERT INTO dataset (dataset_id, section_id, title, sort_order) VALUES (%s, 'other', 'T', 0)",
                        (DATASET,))
            cur.execute("INSERT INTO metric VALUES (%s, 'wet_area', 'Wet area', 'km2', 'sum')", (DATASET,))
            series_id = common.ensure_series(cur, [(DATASET, "wet_area", ())])[(DATASET, "wet_area", ())]
            for boundary_set, rows in ((BASINS, [(keys["1"], 10.0), (keys["2"], 20.0)]), (GROUPS, [(group, 30.0)])):
                load_id = common.begin_load(cur, DATASET, boundary_set, None, "fixture")
                cur.executemany("INSERT INTO observation VALUES (%s, %s, %s)", [(series_id, k, v) for k, v in rows])
                common.finish_load(cur, load_id)
        svc = Service(OneConnection(conn), catalogue(), countries={"KEN", "UGA"})
        svc.refresh()
        yield svc
        conn.rollback()


def code(call):
    with pytest.raises(ApiError) as error:
        call()
    return error.value.status, error.value.code


def test_values_carry_share_and_whole_area_values(service):
    ken = service.values("KEN", "test", "basins", 5, None, None, [])
    uga = service.values("uga", "test", "basins", 5, None, None, [])
    assert ken["boundary_set"] == BASINS                                     # the dataset's first set
    assert [(v["area_id"], v["value"]) for v in ken["values"]] == [("1", 10.0), ("2", 20.0)]
    assert [(v["area_id"], v["value"]) for v in uga["values"]] == [("1", 10.0)]   # same value under both
    assert ken["values"][0]["share"] == pytest.approx(0.5, abs=1e-3)
    assert ken["values"][1]["share"] == pytest.approx(4 / 4.005, abs=1e-4)
    assert code(lambda: service.values("KEN", "test", "basins", 1, None, None, [])) == (404, "unknown_level")


def test_areas_scoped_to_the_country(service):
    body = service.area_observations("KEN", "test", "basins", ["2", "1"], BASINS, None, [])
    assert [(a["area_id"], a["level"], a["observations"][0]["value"]) for a in body["areas"]] == [("2", 5, 20.0), ("1", 5, 10.0)]
    assert code(lambda: service.area_observations("UGA", "test", "basins", ["2"], BASINS, None, [])) == (404, "unknown_area")
    group = service.area_observations("UGA", "test", "basins", ["EAST"], GROUPS, None, [])
    assert group["areas"][0]["share"] is None and group["areas"][0]["observations"][0]["value"] == 30.0
    assert code(lambda: service.area_observations("KEN", "test", "basins", ["KEN"], ADMIN, None, [])) == \
        (404, "unknown_boundary_set")                                       # not one of the dataset's sets


def test_overlay_boundaries_are_whole_with_share(service):
    collection = json.loads(service.boundaries("UGA", 5, BASINS))
    [basin] = collection["features"]
    assert set(basin["properties"]) == {"area_id", "name", "parent_area_id", "share"}
    assert basin["properties"]["share"] == pytest.approx(0.5, abs=1e-3)
    xs = [x for x, _ in basin["geometry"]["coordinates"][0][0]]
    assert min(xs) == 32 and max(xs) == 36                                    # not clipped at the 34°E border
    admin = json.loads(service.boundaries("UGA", 1, None))                   # the country's admin set
    assert admin["features"][0]["properties"] == {"area_id": "UGA.1", "name": "All", "parent_area_id": "UGA"}
    assert code(lambda: service.boundaries("UGA", 0, GROUPS)) == (404, "unknown_boundary_set")
    assert code(lambda: service.boundaries("UGA", 6, BASINS)) == (404, "unknown_level")


def test_country_coverage_and_sets(service):
    uga = service.country("UGA")
    assert uga["boundary_set"] == ADMIN and [lv["level"] for lv in uga["levels"]] == [0, 1]
    assert [s["id"] for s in uga["boundary_sets"]] == [ADMIN, BASINS, GROUPS]
    assert uga["datasets"] == [{"id": DATASET, "boundary_sets": [{"id": BASINS, "levels": [5]},
                                                                 {"id": GROUPS, "levels": [0]}], "download": True}]
    sets = {s["id"]: s for s in service.catalogue_body()["boundary_sets"]}
    assert sets[BASINS]["kind"] == "overlay" and sets[BASINS]["levels"] == [{"level": 5, "label": "Sub-basin"}]


def test_overlay_download_has_share(service):
    _, text = service.download("KEN", "test", "basins", 5, None)
    rows = list(csv.DictReader(io.StringIO(text)))
    assert list(rows[0]) == ["country_iso3", "boundary_set", "level", "area_id", "area_name", "share",
                             "metric", "unit", "value"]
    assert [(r["area_id"], r["value"]) for r in rows] == [("1", "10.0"), ("2", "20.0")]
