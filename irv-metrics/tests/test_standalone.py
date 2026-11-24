"""A deployment that is not GRI: the C2 example `standalone-jam` seeds an empty cube and serves /v1.

The example config (`datasets/examples/standalone-jam/`: one section, one
dataset, a range dimension, its own admin set and country) is loaded with the
CSV loader into a freshly created `metrics` schema, from small generated files
(three GeoPackages and the ministry's CSV). The catalogue is then read back from
the database and through the /v1 service.

The schema is dropped and recreated inside one transaction that is rolled back,
so an existing cube is left as it was. If another session holds locks on the
cube, the test is skipped after a few seconds rather than waiting.

Needs a PostGIS database (DATABASE_URL); skipped without one.
"""
from __future__ import annotations

from contextlib import contextmanager
import csv
import os
from pathlib import Path
import struct

import pytest

from irv_metrics import config
from irv_metrics.config import load_catalogue

pytestmark = pytest.mark.skipif(not os.environ.get("DATABASE_URL"), reason="needs DATABASE_URL (PostGIS)")

EXAMPLE = config.PACKAGE.parent / "datasets" / "examples" / "standalone-jam"
DATASET, SET = "agriculture/crop-production", "jam-statin-2011"


def box(x0, y0, x1, y1) -> bytes:
    """A polygon as WKB."""
    ring = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    return struct.pack("<BIII", 1, 3, 1, len(ring)) + b"".join(struct.pack("<dd", *p) for p in ring)


def write_gpkg(path: Path, rows: list[dict], geometries: list[bytes]):
    import numpy
    import pyogrio.raw

    fields = list(rows[0])
    data = [numpy.array([r[f] for r in rows], dtype=object) for f in fields]
    pyogrio.raw.write(path, numpy.array(geometries, dtype=object), data, fields, driver="GPKG",
                      geometry_type="Polygon", crs="EPSG:4326")


def write_data(root: Path):
    """The files the example's patterns name: country, parishes, enumeration districts, the production table."""
    (root / "boundaries").mkdir(parents=True)
    write_gpkg(root / "boundaries/jam_country.gpkg", [{"iso3": "JAM", "name": "Jamaica"}], [box(0, 0, 4, 2)])
    write_gpkg(root / "boundaries/jam_parishes.gpkg",
               [{"PARISH_CODE": "01", "PARISH": "Kingston"}, {"PARISH_CODE": "02", "PARISH": "St. Andrew"}],
               [box(0, 0, 2, 2), box(2, 0, 4, 2)])
    write_gpkg(root / "boundaries/jam_eds.gpkg",
               [{"ED_ID": "01001", "PARISH_CODE": "01"}, {"ED_ID": "01002", "PARISH_CODE": "01"},
                {"ED_ID": "02001", "PARISH_CODE": "02"}],
               [box(0, 0, 1, 2), box(1, 0, 2, 2), box(2, 0, 4, 2)])
    (root / "agriculture").mkdir()
    with open(root / "agriculture/crop_production_by_parish.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["parish_code", "parish_name", "indicator", "commodity", "crop_year", "amount", "source"])
        for code, name, production, area in (("01", "Kingston", 100.0, 4.0), ("02", "St. Andrew", 250.0, "")):
            for year in ("2023", "2024"):
                writer.writerow([code, name, "production", "all", year, production, "MoA"])
                writer.writerow([code, name, "production", "banana", year, production / 2, "MoA"])
                writer.writerow([code, name, "harvested_area", "all", year, area, "MoA"])


class OneConnection:
    """A pool that hands out one connection (the test's transaction)."""

    def __init__(self, conn):
        self.conn = conn

    @contextmanager
    def connection(self):
        yield self.conn


@pytest.fixture(scope="module")
def cube(tmp_path_factory):
    psycopg = pytest.importorskip("psycopg")
    from irv_metrics.loaders import common
    from irv_metrics.loaders.csv import Loader

    data = tmp_path_factory.mktemp("jam")
    write_data(data)
    catalogue = load_catalogue(EXAMPLE)
    with psycopg.connect(common.dsn()) as conn:
        conn.execute("SET LOCAL lock_timeout = '5s'")
        loader = Loader(conn, catalogue, data, commit=False)
        # Rolled back below, so the reset is allowed where Alembic manages the schema (GRI).
        allowed = os.environ.get(common.ALLOW_RESET_ENV)
        os.environ[common.ALLOW_RESET_ENV] = "1"
        try:
            ok = loader.run(reset=True)                  # drops and recreates the schema, in this transaction
        except psycopg.errors.LockNotAvailable:
            conn.rollback()
            pytest.skip("another session holds locks on the cube")
        finally:
            if allowed is None:
                os.environ.pop(common.ALLOW_RESET_ENV, None)
            else:
                os.environ[common.ALLOW_RESET_ENV] = allowed
        yield conn, catalogue, loader, ok
        conn.rollback()


def test_standalone_config_seeds_and_loads(cube):
    conn, catalogue, loader, ok = cube
    assert ok and not loader.failed, loader.report
    assert loader.report["JAM"]["boundary_set"] == SET
    assert loader.report["JAM"]["areas"] == {0: 1, 1: 2, 2: 3}
    stats = loader.report["JAM"]["datasets"][DATASET]
    assert (stats["placed"], stats["empty"]) == (12, 2)


def test_catalogue_in_the_cube_is_the_config(cube):
    conn, catalogue, *_ = cube
    q = lambda sql: conn.execute(sql).fetchall()  # noqa: E731
    assert q("SELECT section_id, title, sort_order FROM section") == [("agriculture", "Agriculture", 0)]
    assert sorted(q("SELECT unit_id, label, quantity FROM unit")) == [("hectares", "ha", "area"), ("tonnes", "t", "mass")]
    assert q("SELECT dataset_id, section_id, title, source_note FROM dataset") == [
        (DATASET, "agriculture", "Crop production", "Ministry of Agriculture (not open data)")]
    assert q("SELECT metric_id, unit_id, aggregation FROM metric ORDER BY 1") == [
        ("harvested_area", "hectares", "sum"), ("production", "tonnes", "sum")]
    years = q("SELECT member_id, numeric_value, sort_order FROM dimension_member WHERE dimension_id = 'year'"
              " ORDER BY sort_order")
    assert years == [(str(y), float(y), i) for i, y in enumerate(range(2015, 2025))]
    assert q("SELECT is_ordered FROM dimension WHERE dimension_id = 'year'") == [(True,)]
    assert q("SELECT dimension_id, default_member FROM dataset_dimension ORDER BY 1") == [("crop", "all"), ("year", "2024")]
    assert q("SELECT count(*) FROM dataset_dimension_member WHERE dimension_id = 'year'") == [(10,)]
    assert q("SELECT boundary_set, kind, attribution FROM boundary_set") == [
        (SET, "admin", "Statistical Institute of Jamaica")]
    assert q("SELECT level, label FROM boundary_set_level ORDER BY 1") == [
        (0, "Country"), (1, "Parish"), (2, "Enumeration district")]
    assert q("SELECT iso3, name, boundary_set FROM country") == [("JAM", "Jamaica", SET)]
    assert q("SELECT a.area_id, p.area_id FROM area a LEFT JOIN area p ON p.area_key = a.parent_area_key"
             " ORDER BY a.level, a.area_id") == [
        ("JAM", None), ("01", "JAM"), ("02", "JAM"), ("01001", "01"), ("01002", "01"), ("02001", "02")]


def test_v1_service_reads_the_standalone_cube(cube):
    from irv_metrics.api import Service

    conn, catalogue, *_ = cube
    service = Service(OneConnection(conn), catalogue)
    service.refresh()
    body = service.catalogue_body()
    assert [s["id"] for s in body["sections"]] == ["agriculture"]
    assert [d["id"] for d in body["datasets"]] == [DATASET]
    assert [m["id"] for m in body["datasets"][0]["metrics"]] == ["production", "harvested_area"]
    assert body["boundary_sets"][0]["levels"] == [
        {"level": 0, "label": "Country"}, {"level": 1, "label": "Parish"}, {"level": 2, "label": "Enumeration district"}]

    (country,) = service.list_countries()
    assert (country["iso3"], country["name"], country["boundary_set"]) == ("JAM", "Jamaica", SET)
    assert [(lv["level"], lv["label"], lv["areas"]) for lv in country["levels"]] == [
        (0, "Country", 1), (1, "Parish", 2), (2, "Enumeration district", 3)]
    assert service.country("jam")["datasets"] == [
        {"id": DATASET, "boundary_sets": [{"id": SET, "levels": [1]}], "download": True}]

    values = service.values("JAM", "agriculture", "crop-production", 1, None, "production",
                            [("crop", "banana"), ("year", "2024")])
    assert values["unit"] == "tonnes" and values["dims"] == {"crop": "banana", "year": "2024"}
    assert [(v["area_id"], v["name"], v["value"]) for v in values["values"]] == [
        ("01", "Kingston", 50.0), ("02", "St. Andrew", 125.0)]
    observations = service.area_observations("JAM", "agriculture", "crop-production", ["02"], None,
                                             "harvested_area", [("year", "2023")])
    assert [(o["metric"], o["dims"], o["value"]) for o in observations["areas"][0]["observations"]] == [
        ("harvested_area", {"crop": "all", "year": "2023"}, None)]
    name, text = service.download("JAM", "agriculture", "crop-production", 1, None)
    assert name == f"JAM_{SET}_level1_agriculture_crop-production.csv"
    assert len(text.splitlines()) == 1 + 12
