"""The implemented /v1 API has exactly the frozen contract C1 (`irv-metrics/openapi.json`).

`irv_metrics` keeps its own copy of the contract's models (`irv_metrics/api/models.py`),
because GRI vendors the package without pawehe's `backend/contract/`. These tests build
the OpenAPI schema from the implementation and compare it with the committed copy of
the contract (synced from pawehe, see README.md). pawehe's `contract.build --check`
(openapi.json matches contract/v1.py) runs in pawehe only.

The comparison is exact, so it needs the toolchain the contract was generated with:
Python >= 3.13 (HTTP 422 is "Unprocessable Content"), FastAPI >= 0.129 (`ctx` and
`input` in `ValidationError`) and pydantic >= 2.11 (`additionalProperties: true`).

Normalisations, and nothing else:
  1. `info` (title, version, description) is the app's, not the contract's: the bare
     router is mounted on an app with the contract's `info`, and the dev app's `info`
     is not compared.
  2. The dev app also serves non-/v1 paths (legacy routes, /health); only `/v1` paths
     are compared there.
"""
import json
from pathlib import Path

from fastapi import FastAPI

from irv_metrics.api import Service, create_router

BACKEND = Path(__file__).resolve().parents[1]
CONTRACT = json.loads((BACKEND / "openapi.json").read_text())


def roundtrip(schema: dict) -> dict:
    return json.loads(json.dumps(schema))


def first_difference(expected, actual, path="$"):
    if type(expected) is not type(actual):
        return f"{path}: {type(expected).__name__} != {type(actual).__name__}"
    if isinstance(expected, dict):
        for key in sorted(set(expected) | set(actual)):
            if key not in actual or key not in expected:
                return f"{path}.{key}: only in {'contract' if key in expected else 'implementation'}"
            found = first_difference(expected[key], actual[key], f"{path}.{key}")
            if found:
                return found
        return None
    if isinstance(expected, list):
        if len(expected) != len(actual):
            return f"{path}: {len(expected)} items != {len(actual)}"
        for i, (a, b) in enumerate(zip(expected, actual)):
            found = first_difference(a, b, f"{path}[{i}]")
            if found:
                return found
        return None
    return None if expected == actual else f"{path}: {expected!r} != {actual!r}"


def router_schema() -> dict:
    app = FastAPI(**{k: CONTRACT["info"][k] for k in ("title", "version", "description") if k in CONTRACT["info"]})
    app.include_router(create_router(Service(pool=None)))
    return roundtrip(app.openapi())


def test_router_schema_equals_contract():
    schema = router_schema()
    assert first_difference(CONTRACT, schema) is None
    assert schema == CONTRACT


def test_operations_are_the_contracts():
    schema = router_schema()
    operations = {(path, method): op["operationId"] for path, item in schema["paths"].items() for method, op in item.items()}
    assert sorted(operations.values()) == sorted(
        ["listCountries", "getCountry", "getCatalogue", "getBoundaries", "getValues", "getAreaObservations",
         "getCurve", "downloadDataset"])


def test_dev_app_serves_the_contract_under_v1():
    from irv_metrics.adapters.fastapi import app     # the pool is created closed; nothing connects

    schema = roundtrip(app.openapi())
    v1 = {path: item for path, item in schema["paths"].items() if path.startswith("/v1/") or path == "/v1"}
    assert first_difference(CONTRACT["paths"], v1, "$.paths") is None
    assert first_difference(CONTRACT["components"], schema["components"], "$.components") is None


def test_package_does_not_import_contract():
    offenders = [str(p.relative_to(BACKEND)) for p in (BACKEND / "irv_metrics").rglob("*.py")
                 if any(line.lstrip().startswith(("from contract", "import contract"))
                        for line in p.read_text().splitlines())]
    assert offenders == []
