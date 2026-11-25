"""The CSV loader on a fixture: discovery and boundary reading always; the load itself with a database.

The database tests need a PostGIS database (DATABASE_URL) and are skipped
without one. They run in one transaction that is rolled back, so an existing
cube is left untouched (the schema is created inside it if missing). The
fixture is described in tests/fixtures/csv/README.md.
"""
from __future__ import annotations

import os
from pathlib import Path

import pytest

from irv_metrics.config import load_catalogue
from irv_metrics.loaders.csv import COMPUTED, Bands, Loader, computed_problems, discover, read_boundaries

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "csv"
CONFIG, DATA = FIXTURE / "config", FIXTURE / "data"
SET = "test-admin"
FLOOD, CURVES = "hazard/river-flooding", "adaptation-outcomes/concentration-curves"


@pytest.fixture(scope="module")
def catalogue():
    return load_catalogue(CONFIG)


def test_discovery_skips_blank_admin0_stubs_and_reports_completeness_gaps(catalogue):
    found = discover(catalogue, DATA)
    assert found.countries() == ["ZZA", "ZZB"]
    assert [s.rel for s in found.stubs] == ["ZZB/hazard/ZZB_adm0_flood.csv"]
    assert found.gaps == [f"ZZB {FLOOD}: admin levels [1] without usable admin 0"]
    assert sorted((f.dataset, f.iso3, f.level) for f in found.of("ZZA")) == [
        (CURVES, "ZZA", None), (FLOOD, "ZZA", 0), (FLOOD, "ZZA", 1)]


def test_boundaries_are_read_for_one_country_and_one_level_type(catalogue):
    areas = read_boundaries(catalogue, SET, DATA, "ZZA")
    assert [(a.level, a.area_id, a.name) for a in areas] == [
        (0, "ZZA", "Zedland"), (1, "ZZA.1", "West"), (1, "ZZA.2", "ZZA.2")]   # unnamed: the id; ADM0 in admin1 dropped
    assert all(isinstance(a.geometry, bytes) for a in areas)


WIND = "hazard/tropical-cyclone-wind"


def test_wind_bands_are_declared_to_fit_the_gri_config():
    gri = load_catalogue()
    assert computed_problems(gri) == []
    assert list(COMPUTED) == [WIND]


def test_bands_are_differences_of_consecutive_thresholds():
    wind = load_catalogue().datasets[WIND]
    (rule,) = COMPUTED[WIND]

    def key(metric, dim, member, rp="10"):
        return metric, tuple(sorted({"return_period_years": rp, dim: member}.items()))

    thresholds = ["tropical_storm_plus", "cat1plus", "cat2plus", "cat3plus", "cat4plus", "cat5plus"]
    area = {1: [10.0, 6.0, 3.5, 1.0, 0.25, 0.0], 2: [5.0, None, 1.0, 1.0, 0.0, 0.0]}   # area_key -> by threshold
    values = {key("area_exceeding", "wind_threshold", t): {k: v[i] for k, v in area.items()}
              for i, t in enumerate(thresholds)}
    values[key("area_exceeding", "wind_threshold", "cat5plus")].pop(2)                  # no top threshold on area 2
    values[key("area_exceeding_share", "wind_threshold", "tropical_storm_plus", "100")] = {1: 50.0}
    out = rule.compute(wind, values)
    assert out[key("area_in_band", "wind_band", "tropical_storm")] == {1: 4.0, 2: None}
    assert out[key("area_in_band", "wind_band", "cat1")] == {1: 2.5, 2: None}
    assert out[key("area_in_band", "wind_band", "cat3")] == {1: 0.75, 2: 1.0}
    assert out[key("area_in_band", "wind_band", "cat4")] == {1: 0.25}                   # area 2: cat5plus missing
    assert out[key("area_in_band", "wind_band", "cat5")] == {1: 0.0}                    # the top band: the threshold
    # A band whose upper threshold is not published at all is not computed.
    assert key("area_in_band_share", "wind_band", "tropical_storm", "100") not in out


def test_band_declarations_are_checked():
    wind = load_catalogue().datasets[WIND]
    wrong = Bands({"area_in_band": "area_exceeding_share"}, threshold="wind_threshold", band="wind_band")
    assert wrong.check(wind) == [f"{WIND}: area_in_band and area_exceeding_share differ in unit"]
    assert Bands({}, "wind_threshold", "sector").check(wind) == [f"{WIND}: needs dimensions wind_threshold and sector"]


# --- with a database -----------------------------------------------------------

needs_db = pytest.mark.skipif(not os.environ.get("DATABASE_URL"), reason="needs DATABASE_URL (PostGIS)")


def snapshot(cur) -> dict:
    cur.execute("""SELECT s.dataset_id, s.metric_id, s.dims, a.area_id, a.level, o.value
                   FROM observation o JOIN series s USING (series_id) JOIN area a USING (area_key)
                   WHERE a.boundary_set = %s""", (SET,))
    values = {(d, m, tuple(sorted(dims.items())), area): (level, v) for d, m, dims, area, level, v in cur.fetchall()}
    cur.execute("""SELECT s.dims->>'scenario', p.seq, p.x, p.y FROM curve_point p JOIN series s USING (series_id)
                   JOIN area a USING (area_key) WHERE a.boundary_set = %s ORDER BY 1, 2""", (SET,))
    curves = cur.fetchall()
    cur.execute("""SELECT a.area_id, a.name, p.area_id, ac.iso3 FROM area a LEFT JOIN area p ON p.area_key = a.parent_area_key
                   JOIN area_country ac ON ac.area_key = a.area_key WHERE a.boundary_set = %s ORDER BY 1""", (SET,))
    areas = cur.fetchall()
    cur.execute("SELECT dataset_id, iso3, finished_at IS NOT NULL FROM load_record WHERE boundary_set = %s ORDER BY 1",
                (SET,))
    records = cur.fetchall()
    cur.execute("SELECT name, boundary_set FROM country WHERE iso3 = 'ZZA'")
    return {"values": values, "curves": curves, "areas": areas, "records": records, "country": cur.fetchone()}


@pytest.fixture(scope="module")
def loaded(catalogue):
    psycopg = pytest.importorskip("psycopg")
    from irv_metrics.loaders import common

    with psycopg.connect(common.dsn()) as conn:
        loader = Loader(conn, catalogue, DATA, commit=False)
        ok = loader.run(["ZZA"])
        with conn.cursor() as cur:
            first = snapshot(cur)
        yield conn, loader, ok, first
        conn.rollback()


@needs_db
def test_fixture_load(loaded):
    _, loader, ok, cube = loaded
    assert ok and not loader.failed
    rp = lambda n: (("return_period_years", n),)  # noqa: E731
    assert cube["values"] == {
        (FLOOD, "flooded_area", rp("10"), "ZZA"): (0, 12.5),       # adm_id "Zedland" is read as ZZA
        (FLOOD, "flooded_area", rp("100"), "ZZA"): (0, 30.25),
        (FLOOD, "flooded_area", rp("10"), "ZZA.1"): (1, 5.0),
        (FLOOD, "flooded_area", rp("10"), "ZZA.2"): (1, 7.5),
        (FLOOD, "flooded_area", rp("100"), "ZZA.1"): (1, None),    # published without a value: NULL
        (FLOOD, "flooded_area", rp("100"), "ZZA.2"): (1, 20.0),
    }
    assert cube["curves"] == [("baseline_protected", 0, 0.0, 0.0), ("baseline_protected", 1, 0.5, 0.7),
                              ("baseline_protected", 2, 1.0, 1.0), ("dry_proofing", 0, 0.0, 0.0),
                              ("dry_proofing", 1, 0.5, 0.6), ("dry_proofing", 2, 1.0, 1.0)]
    assert cube["areas"] == [("ZZA", "Zedland", None, "ZZA"), ("ZZA.1", "West", "ZZA", "ZZA"),
                             ("ZZA.2", "ZZA.2", "ZZA", "ZZA")]
    assert cube["records"] == [(CURVES, "ZZA", True), (FLOOD, "ZZA", True)]
    assert cube["country"] == ("Zedland", SET)
    assert loader.report["ZZA"]["datasets"][FLOOD]["empty"] == 1


@needs_db
def test_skip_and_replace_leave_the_cube_unchanged(loaded, catalogue):
    conn, _, _, first = loaded
    skip = Loader(conn, catalogue, DATA, if_loaded="skip", commit=False)
    assert skip.run(["ZZA"]) and skip.report["ZZA"].get("skipped")
    with conn.cursor() as cur:
        assert snapshot(cur) == first
    replace = Loader(conn, catalogue, DATA, if_loaded="replace", commit=False)
    assert replace.run(["ZZA"]) and not replace.report["ZZA"].get("skipped")
    with conn.cursor() as cur:
        assert snapshot(cur) == first


@needs_db
def test_completeness_gap_fails_before_writing(loaded, catalogue):
    conn = loaded[0]
    loader = Loader(conn, catalogue, DATA, commit=False)
    assert not loader.run(["ZZB"])
    with conn.cursor() as cur:
        cur.execute("SELECT count(*) FROM area_country WHERE iso3 = 'ZZB'")
        assert cur.fetchone()[0] == 0
