"""The /v1 routes (contract C1) as a FastAPI router.

    from irv_metrics.api import Service, create_router
    service = Service(pool)                       # pool: psycopg_pool.ConnectionPool, search_path metrics
    app.include_router(create_router(service), prefix="/api/national-metrics")

The route declarations (paths, parameters, operation ids, docstrings, models)
are those of `backend/contract/v1.py`; `tests/test_contract_conformance.py`
checks that the generated OpenAPI matches `backend/openapi.json`.

Everything a host would otherwise configure on its app is done by the route
class, so the router behaves the same wherever it is mounted:
  * errors, including request validation errors, are `{"detail": {"code", "message"}}`;
  * pool exhaustion and statement timeouts are 503 `busy`;
  * every 200 carries `ETag` (the catalogue version) and `Cache-Control`, and a
    matching `If-None-Match` is answered with 304.
"""
from __future__ import annotations

from typing import Annotated, Callable

from fastapi import APIRouter, Query, Request
from fastapi.concurrency import run_in_threadpool
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse, Response
from fastapi.routing import APIRoute
from fastapi import Path
import psycopg
from psycopg_pool import PoolTimeout

from .models import (DIMS, FILTERED, NOT_FOUND, AreaIds, AreaObservations, BoundarySet, Catalogue, Country,
                     CountryDetail, Curve, Error, DatasetName, Iso3, Level, Metric, Section, SeriesValues)
from .service import ApiError, Service

HOUR = 3600
DATA_MAX_AGE = 300
# Pool exhausted; statement timeout (QueryCanceled), too many connections or the database gone (OperationalError).
BUSY = (PoolTimeout, psycopg.OperationalError)


def max_age(seconds: int):
    """Mark a route's Cache-Control max-age."""
    def mark(endpoint):
        endpoint.max_age = seconds
        return endpoint
    return mark


def error(status: int, code: str, message: str) -> JSONResponse:
    return JSONResponse({"detail": {"code": code, "message": message}}, status_code=status,
                        headers={"Cache-Control": "no-store"})


def validation_error(exc: RequestValidationError) -> JSONResponse:
    problems = list(exc.errors())
    for problem in problems:
        loc = list(problem.get("loc", ()))
        if loc[:2] == ["path", "iso3"]:
            return error(404, "unknown_country", f"Not a country code: {problem.get('input')!r}")
        if loc[:2] == ["path", "level"]:
            return error(404, "unknown_level", f"Not a level: {problem.get('input')!r}")
    message = "; ".join(f"{'.'.join(str(p) for p in problem.get('loc', ())[1:]) or 'request'}: {problem.get('msg')}"
                        for problem in problems)
    return error(422, "invalid_filters", message)


def etag_matches(header: str | None, etag: str) -> bool:
    if not header:
        return False
    if header.strip() == "*":
        return True
    return any(tag.strip().removeprefix("W/") == etag for tag in header.split(","))


def route_class(service: Service) -> type[APIRoute]:
    class V1Route(APIRoute):
        def get_route_handler(self) -> Callable:
            inner = super().get_route_handler()
            seconds = getattr(self.endpoint, "max_age", DATA_MAX_AGE)

            async def handler(request: Request) -> Response:
                try:
                    await run_in_threadpool(service.ensure_fresh)
                    etag = f'"{service.current().version}"'
                    response = await inner(request)
                except ApiError as exc:
                    return error(exc.status, exc.code, exc.message)
                except RequestValidationError as exc:
                    return validation_error(exc)
                except BUSY as exc:
                    return error(503, "busy", f"The database is busy; retry shortly ({type(exc).__name__})")
                if response.status_code != 200:
                    return response
                headers = {"ETag": etag, "Cache-Control": f"public, max-age={seconds}"}
                if etag_matches(request.headers.get("if-none-match"), etag):
                    return Response(status_code=304, headers=headers)
                response.headers.update(headers)
                return response

            return handler

    return V1Route


def create_router(service: Service) -> APIRouter:
    """The /v1 router over a `Service`. The host refreshes the service once its pool is open
    (`service.refresh()`); after that each request checks for a newer load."""
    router = APIRouter(prefix="/v1", route_class=route_class(service))

    # ------------------------------------------------------------------ catalogue

    @router.get("/countries", response_model=list[Country], operation_id="listCountries", tags=["catalogue"])
    @max_age(HOUR)
    def list_countries():
        """Countries with data. A single-country deployment returns one."""
        return service.list_countries()

    @router.get("/countries/{iso3}", response_model=CountryDetail, operation_id="getCountry",
                responses=NOT_FOUND, tags=["catalogue"])
    @max_age(HOUR)
    def get_country(iso3: Iso3):
        """One country and its coverage: which datasets have data, at which levels."""
        return service.country(iso3)

    @router.get("/datasets", response_model=Catalogue, operation_id="getCatalogue", tags=["catalogue"])
    @max_age(HOUR)
    def get_catalogue():
        """Everything needed to build controls and labels: sections, units, dimensions with
        their members, and each dataset's metrics, dimensions and display hints.

        The same for every country. Use `getCountry` for what a country actually has.
        """
        return service.catalogue_body()

    # ------------------------------------------------------------------ geometry

    @router.get("/countries/{iso3}/boundaries/{level}", operation_id="getBoundaries", tags=["geometry"],
                responses={200: {"content": {"application/geo+json": {"schema": {"type": "object"}}},
                                 "description": "GeoJSON FeatureCollection. Each feature has `id` = `area_id` and "
                                                "properties `area_id`, `name`, `parent_area_id`, and `share` in "
                                                "overlay sets. Overlay areas are not clipped to the border."},
                           404: {"model": Error}},
                response_class=Response)
    @max_age(HOUR)
    def get_boundaries(iso3: Iso3, level: Annotated[int, Path(ge=0)], boundary_set: BoundarySet = None):
        """Geometry of the areas at one level related to the country, WGS84, 6 decimal places.

        Group sets have no geometry: 404 `unknown_boundary_set`.
        """
        return Response(service.boundaries(iso3, level, boundary_set), media_type="application/geo+json")

    # ------------------------------------------------------------------ data

    @router.get("/countries/{iso3}/datasets/{section}/{dataset}/values", response_model=SeriesValues,
                operation_id="getValues", responses=FILTERED, tags=["data"], openapi_extra={"parameters": [DIMS]})
    def get_values(request: Request, iso3: Iso3, section: Section, dataset: DatasetName, level: Level,
                   boundary_set: BoundarySet = None, metric: Metric = None):
        """One series across every area at a level: the data behind a choropleth.

        The filters must be the complete dimension assignment of one stored series. The
        server fills in no defaults; the catalogue gives each dimension's default.

        422 `invalid_filters` when they are not: an unknown dimension or member, or an
        assignment that is no series (usually a dimension left out). The message says which.
        """
        return service.values(iso3, section, dataset, level, boundary_set, metric,
                              request.query_params.multi_items())

    @router.get("/countries/{iso3}/datasets/{section}/{dataset}/areas", response_model=AreaObservations,
                operation_id="getAreaObservations", responses=FILTERED, tags=["data"],
                openapi_extra={"parameters": [DIMS]})
    def get_area_observations(request: Request, iso3: Iso3, section: Section, dataset: DatasetName,
                              area_id: AreaIds, boundary_set: BoundarySet = None, metric: Metric = None):
        """Every stored value for some areas of one set that match the filters.

        Filters are partial: a dimension left out is left free. This one route serves what
        the upstream API split into "profile" (all members of a dimension) and "ladder"
        (all return periods). Several areas in one call serve comparisons, e.g. a county and
        its country. Areas in another set (a country and its income group) take a second call.

        404 `unknown_area` when any `area_id` is not in the set or not related to the country.

        A filter on a dimension that a metric is not keyed by does not exclude that metric.
        Derived metrics are returned only when `metric` names one. A valid filter that
        matches nothing returns an empty list.
        """
        return service.area_observations(iso3, section, dataset, area_id, boundary_set, metric,
                                         request.query_params.multi_items())

    @router.get("/countries/{iso3}/datasets/{section}/{dataset}/curve", response_model=Curve,
                operation_id="getCurve", responses=FILTERED, tags=["data"], openapi_extra={"parameters": [DIMS]})
    def get_curve(request: Request, iso3: Iso3, section: Section, dataset: DatasetName,
                  boundary_set: BoundarySet = None, metric: Metric = None,
                  area_id: Annotated[str | None, Query(description="Defaults to the country, in an admin set.")] = None):
        """The points of one curve, for datasets of kind `curve`. Filters are complete, as for `getValues`."""
        return service.curve(iso3, section, dataset, boundary_set, metric, area_id,
                             request.query_params.multi_items())

    @router.get("/countries/{iso3}/datasets/{section}/{dataset}/download", operation_id="downloadDataset",
                tags=["data"], response_class=Response,
                responses={200: {"content": {"text/csv": {"schema": {"type": "string"}}},
                                 "description": "One row per area and series: key columns, one column per dimension, "
                                                "`metric`, `unit`, `value`."},
                           404: {"model": Error}})
    def download_dataset(iso3: Iso3, section: Section, dataset: DatasetName, level: Level,
                         boundary_set: BoundarySet = None):
        """The dataset for one country and level as CSV, generated from the cube.

        404 `downloads_disabled` where the deployment turns downloads off.
        """
        filename, text = service.download(iso3, section, dataset, level, boundary_set)
        return Response(text, media_type="text/csv",
                        headers={"Content-Disposition": f'attachment; filename="{filename}"'})

    return router

