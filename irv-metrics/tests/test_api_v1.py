"""The /v1 routes over the loaded cube (KEN from the scrape): semantics, errors and caching.

Needs DATABASE_URL with KEN loaded (`python -m irv_metrics.loaders.scrape ../api-data/KEN`),
or IRV_METRICS_BASE_URL; skipped otherwise. The busy (503) tests need no database.
"""
import csv
import io
import json
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import psycopg
from psycopg_pool import PoolTimeout
import pytest
from fastapi import FastAPI

from irv_metrics.api import Service, create_router
from irv_metrics.config import load_catalogue

CATALOGUE = load_catalogue()


def get(base, path, params=(), headers=None):
    url = base + path + ("?" + urlencode(list(params)) if params else "")
    try:
        with urlopen(Request(url, headers=headers or {}), timeout=60) as response:
            return response.status, response.read(), response.headers
    except HTTPError as error:
        return error.code, error.read(), error.headers


def ok(base, path, params=()):
    status, body, headers = get(base, path, params)
    assert status == 200, body
    return json.loads(body) if "json" in headers.get("content-type", "") else body


def problem(base, path, params=()):
    status, body, headers = get(base, path, params)
    detail = json.loads(body)
    assert set(detail) == {"detail"} and set(detail["detail"]) == {"code", "message"}, detail
    assert headers["Cache-Control"] == "no-store"
    return status, detail["detail"]["code"]


@pytest.fixture(scope="module")
def base(api_base):
    status, body, _ = get(api_base, "/v1/countries/KEN")
    if status != 200:
        pytest.skip("KEN is not loaded")
    return api_base


KEN = "/v1/countries/KEN"


# ---------------------------------------------------------------- catalogue

def test_countries_and_country(base):
    countries = {c["iso3"]: c for c in ok(base, "/v1/countries")}
    ken = countries["KEN"]
    west, south, east, north = ken["bbox"]
    assert 33 < west < east < 43 and -5 < south < north < 6
    assert ken["boundary_set"] == "cgaz-6.0.0"
    assert [lv["level"] for lv in ken["levels"]] == [0, 1, 2] and ken["levels"][0]["areas"] == 1
    detail = ok(base, "/v1/countries/ken")
    assert detail["iso3"] == "KEN" and detail["boundary_sets"][0]["id"] == "cgaz-6.0.0"
    coverage = {d["id"]: d for d in detail["datasets"]}
    assert coverage["hazard/river-flooding"]["boundary_sets"] == [{"id": "cgaz-6.0.0", "levels": [0, 1, 2]}]
    assert coverage["adaptation-outcomes/concentration-curves"]["boundary_sets"][0]["levels"] == [0]
    order = list(CATALOGUE.datasets)
    assert [order.index(d) for d in coverage] == sorted(order.index(d) for d in coverage)


def test_catalogue_is_the_config(base):
    body = ok(base, "/v1/datasets")
    assert [d["id"] for d in body["datasets"]] == list(CATALOGUE.datasets)
    population = next(d for d in body["datasets"] if d["id"] == "exposure/population")
    share = next(m for m in population["metrics"] if m["id"] == "population_share")
    assert share["derived"] is True and share["dimensions"] == ["demographic_group"] and share["unit"] == "percent"
    assert population["display"]["card"] == "population"
    curves = next(d for d in body["datasets"] if d["id"] == "adaptation-outcomes/concentration-curves")
    assert curves["kind"] == "curve" and curves["x"]["id"] == "cumulative_population_share"
    assert [s["id"] for s in body["boundary_sets"]] == list(CATALOGUE.boundary_sets)
    assert len(body["version"]) == 16


def test_etag_cache_control_and_304(base):
    status, body, headers = get(base, "/v1/datasets")
    version = json.loads(body)["version"]
    assert headers["ETag"] == f'"{version}"'
    assert headers["Cache-Control"] == "public, max-age=3600"
    status, body, headers = get(base, "/v1/datasets", headers={"If-None-Match": f'"{version}"'})
    assert status == 304 and body == b"" and headers["ETag"] == f'"{version}"'
    assert get(base, "/v1/datasets", headers={"If-None-Match": '"other"'})[0] == 200
    status, _, headers = get(base, f"{KEN}/datasets/exposure/population/values",
                             [("level", 1), ("demographic_group", "total")])
    assert status == 200 and headers["ETag"] == f'"{version}"' and headers["Cache-Control"] == "public, max-age=300"
    assert get(base, f"{KEN}/boundaries/1")[2]["Cache-Control"] == "public, max-age=3600"


# ---------------------------------------------------------------- values

def test_values_one_series(base):
    body = ok(base, f"{KEN}/datasets/risk/capital-stock/values",
              [("level", 1), ("metric", "loss"), ("sector", "total"), ("return_period_years", "100")])
    assert body["dims"] == {"sector": "total", "return_period_years": "100"}
    assert body["unit"] == "usd" and body["boundary_set"] == "cgaz-6.0.0" and body["level"] == 1
    ids = [v["area_id"] for v in body["values"]]
    assert len(ids) == 47 and ids == sorted(ids)
    assert all(v["share"] is None for v in body["values"])


def test_values_ignores_dimensions_the_metric_is_not_keyed_by(base):
    body = ok(base, f"{KEN}/datasets/risk/capital-stock/values",
              [("level", 1), ("metric", "average_annual_loss"), ("sector", "total"), ("return_period_years", "100")])
    assert body["dims"] == {"sector": "total"}


def test_values_derived_share(base):
    path = f"{KEN}/datasets/vulnerability/wealth-distribution/values"
    shares = [ok(base, path, [("level", 1), ("metric", "population_share"), ("wealth_group", q)])
              for q in ("q1", "q2", "q3", "q4", "q5")]
    assert shares[0]["unit"] == "percent" and shares[0]["metric"] == "population_share"
    for i, area in enumerate(shares[0]["values"]):
        assert abs(sum(s["values"][i]["value"] for s in shares) - 100) < 1e-9


@pytest.mark.parametrize("params, status, code", [
    ([("level", 1), ("metric", "loss"), ("sector", "total")], 422, "invalid_filters"),            # rp left out
    ([("level", 1), ("metric", "loss"), ("sector", "everything"), ("return_period_years", "100")], 422, "invalid_filters"),
    ([("level", 1), ("metric", "loss"), ("colour", "red")], 422, "invalid_filters"),              # unknown dimension
    ([("level", 1), ("sector", "total")], 422, "invalid_filters"),                                # metric required
    ([("level", 1), ("metric", "nope"), ("sector", "total")], 422, "invalid_filters"),
    ([("level", 1), ("metric", "loss"), ("sector", "total"), ("sector", "residential"),
      ("return_period_years", "100")], 422, "invalid_filters"),                                   # repeated
    ([("metric", "average_annual_loss"), ("sector", "total")], 422, "invalid_filters"),            # level missing
    ([("level", "x"), ("metric", "average_annual_loss"), ("sector", "total")], 422, "invalid_filters"),
    ([("level", 7), ("metric", "average_annual_loss"), ("sector", "total")], 404, "unknown_level"),
    ([("level", 1), ("set", "hydrobasins"), ("metric", "average_annual_loss"), ("sector", "total")], 404,
     "unknown_boundary_set"),
])
def test_values_errors(base, params, status, code):
    assert problem(base, f"{KEN}/datasets/risk/capital-stock/values", params) == (status, code)


def test_unknown_country_and_dataset(base):
    assert problem(base, "/v1/countries/XXX/datasets/risk/capital-stock/values", [("level", 1)]) == (404, "unknown_country")
    assert problem(base, "/v1/countries/KE") == (404, "unknown_country")
    assert problem(base, f"{KEN}/datasets/risk/nothing/values", [("level", 1)]) == (404, "unknown_dataset")
    assert problem(base, f"{KEN}/datasets/adaptation-outcomes/concentration-curves/values",
                   [("level", 0)]) == (404, "unknown_dataset")


# ---------------------------------------------------------------- areas

def test_areas_two_ids_in_request_order(base):
    county = ok(base, f"{KEN}/datasets/risk/capital-stock/values",
                [("level", 1), ("metric", "average_annual_loss"), ("sector", "total")])["values"][0]["area_id"]
    body = ok(base, f"{KEN}/datasets/risk/capital-stock/areas",
              [("area_id", county), ("area_id", "KEN"), ("metric", "loss"), ("sector", "total")])
    assert [a["area_id"] for a in body["areas"]] == [county, "KEN"]
    assert [a["level"] for a in body["areas"]] == [1, 0]
    rps = [o["dims"]["return_period_years"] for o in body["areas"][1]["observations"]]
    assert rps == [m for m in CATALOGUE.members("risk/capital-stock", "return_period_years") if m in rps]
    assert len(rps) > 3 and all(o["metric"] == "loss" and o["unit"] == "usd" for o in body["areas"][1]["observations"])


def test_areas_order_metric_then_members(base):
    body = ok(base, f"{KEN}/datasets/risk/capital-stock/areas", [("area_id", "KEN")])
    observations = body["areas"][0]["observations"]
    metrics = [o["metric"] for o in observations]
    assert metrics == sorted(metrics, key=list(CATALOGUE.datasets["risk/capital-stock"].metrics).index)
    assert "population_share" not in {o["metric"] for o in
                                      ok(base, f"{KEN}/datasets/exposure/population/areas", [("area_id", "KEN")])
                                      ["areas"][0]["observations"]}


def test_areas_filter_on_unkeyed_dimension_keeps_metric(base):
    body = ok(base, f"{KEN}/datasets/exposure/healthcare-facilities/areas",
              [("area_id", "KEN"), ("demographic_group", "female")])
    found = {(o["metric"], tuple(sorted(o["dims"].items()))) for o in body["areas"][0]["observations"]}
    assert ("count", ()) in found and ("per_100k", (("demographic_group", "female"),)) in found
    assert all(dict(d).get("demographic_group", "female") == "female" for _, d in found)


def test_areas_derived_only_when_asked(base):
    body = ok(base, f"{KEN}/datasets/exposure/population/areas", [("area_id", "KEN"), ("metric", "population_share")])
    observations = body["areas"][0]["observations"]
    assert {o["metric"] for o in observations} == {"population_share"}
    total = next(o for o in observations if o["dims"] == {"demographic_group": "total"})
    assert total["value"] == pytest.approx(100) and total["unit"] == "percent"


def test_areas_errors(base):
    path = f"{KEN}/datasets/risk/capital-stock/areas"
    assert problem(base, path, [("area_id", "KEN"), ("area_id", "MOZ")]) == (404, "unknown_area")
    assert problem(base, path, [("area_id", "nowhere")]) == (404, "unknown_area")
    assert problem(base, path, []) == (422, "invalid_filters")                                      # no area_id
    assert problem(base, path, [("area_id", "KEN")] * 51) == (422, "invalid_filters")               # over 50
    # Valid filters that match nothing: risk/population series name a demographic or a wealth group, not both.
    assert ok(base, f"{KEN}/datasets/risk/population/areas",
              [("area_id", "KEN"), ("demographic_group", "female"), ("wealth_group", "q1")]
              )["areas"][0]["observations"] == []


# ---------------------------------------------------------------- curve, boundaries, download

def test_curve(base):
    path = f"{KEN}/datasets/adaptation-outcomes/concentration-curves/curve"
    body = ok(base, path, [("metric", "cumulative_exposure_share"), ("scenario", "baseline_protected")])
    assert body["area_id"] == "KEN" and body["dims"] == {"scenario": "baseline_protected"}
    assert body["x"]["id"] == "cumulative_population_share" and len(body["points"]) > 10
    assert body["points"][0] == {"x": 0.0, "y": 0.0}
    assert problem(base, path, [("scenario", "baseline_protected"), ("hazard", "river_flood")]) == (422, "invalid_filters")
    assert problem(base, path, [("scenario", "relocation")]) == (422, "invalid_filters")         # threshold left out
    assert problem(base, f"{KEN}/datasets/risk/capital-stock/curve", []) == (404, "unknown_dataset")


def test_boundaries(base):
    status, body, headers = get(base, f"{KEN}/boundaries/2")
    assert status == 200 and headers["content-type"] == "application/geo+json"
    collection = json.loads(body)
    counties = {f["id"] for f in json.loads(get(base, f"{KEN}/boundaries/1")[1])["features"]}
    feature = collection["features"][0]
    assert set(feature["properties"]) == {"area_id", "name", "parent_area_id"}
    assert feature["id"] == feature["properties"]["area_id"]
    assert all(f["properties"]["parent_area_id"] in counties for f in collection["features"])
    assert problem(base, f"{KEN}/boundaries/5") == (404, "unknown_level")
    assert problem(base, f"{KEN}/boundaries/1", [("set", "nope")]) == (404, "unknown_boundary_set")


def test_download(base):
    status, body, headers = get(base, f"{KEN}/datasets/risk/capital-stock/download", [("level", 1)])
    assert status == 200 and headers["content-type"].startswith("text/csv")
    assert "attachment" in headers["Content-Disposition"]
    rows = list(csv.DictReader(io.StringIO(body.decode())))
    assert list(rows[0]) == ["country_iso3", "boundary_set", "level", "area_id", "area_name",
                             "sector", "return_period_years", "metric", "unit", "value"]
    series = {(r["metric"], r["sector"], r["return_period_years"]) for r in rows}
    assert len(rows) == 47 * len(series)
    one = ok(base, f"{KEN}/datasets/risk/capital-stock/values",
             [("level", 1), ("metric", "loss"), ("sector", "total"), ("return_period_years", "100")])["values"][0]
    row = next(r for r in rows if r["area_id"] == one["area_id"] and r["metric"] == "loss"
               and r["sector"] == "total" and r["return_period_years"] == "100")
    assert float(row["value"]) == one["value"]
    assert problem(base, f"{KEN}/datasets/risk/capital-stock/download", [("level", 9)]) == (404, "unknown_level")


# ---------------------------------------------------------------- busy (no database)

class BrokenPool:
    def __init__(self, error):
        self.error = error

    def connection(self):
        raise self.error


@pytest.mark.parametrize("error", [PoolTimeout("pool exhausted"), psycopg.errors.QueryCanceled("statement timeout")])
def test_busy(serve_app, error):
    app = FastAPI()
    app.include_router(create_router(Service(BrokenPool(error))))
    base = serve_app(app)
    assert problem(base, "/v1/datasets") == (503, "busy")
