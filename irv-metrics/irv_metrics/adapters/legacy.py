"""Upstream-compatible routes: the GRI API's paths, parameters and JSON layouts,
rendered from the cube through irv_metrics.specs.

    GET /boundaries/{iso3}/{level}
    GET /metrics/{family}/{layer}/{iso3}/{level}?...
    GET /metrics/{family}/{layer}/{iso3}/profile/{adm_id}?...
    GET /metrics/{family}/{layer}/{iso3}/ladder/{adm_id}?...
    GET /metrics/concentration-curves/{iso3}?...
    GET /downloads/{ISO3}/{ISO3}_adm{1|2}_{card_slug}_metrics.csv

The prototype web app uses this until it is ported to /v1, when we can delete
this module.
"""
from __future__ import annotations

import csv
import io
import re
import time

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import Response

from .. import queries
from ..specs import CURVE_ECHO, CURVES, DATASETS, Dataset, MappingError, fill


class Catalog:
    """Small, rarely changing tables held in memory: series keys, units, areas, member order."""

    version = None
    checked = 0.0

    def __init__(self, pool):
        self.pool = pool

    def refresh(self):
        with self.pool.connection() as conn:
            self.version = queries.load_version(conn)
            self.series = queries.series_keys(conn)
            self.units = queries.metric_units(conn)
            # (boundary_set, area_id) -> (iso3, level, name); the upstream routes only know the CGAZ set.
            self.areas = {(queries.CGAZ, a): v for a, v in queries.admin_areas(conn).items()}
            self.order = queries.member_order(conn)
            self.countries = queries.countries(conn)
        self.boundaries = {}
        self.checked = time.monotonic()

    def ensure_fresh(self):
        """Reload the in-memory catalog when the loader has finished a new load (checked every 5 s)."""
        if time.monotonic() - self.checked < 5:
            return
        self.checked = time.monotonic()
        with self.pool.connection() as conn:
            version = queries.load_version(conn)
        if version != self.version:
            self.refresh()

    def series_id(self, dataset_id: str, metric: str, dims: dict) -> int | None:
        return self.series.get((dataset_id, metric, frozenset(dims.items())))


def fail(status: int, code: str, message: str):
    raise HTTPException(status, {"code": code, "message": message})


def dataset_for(family: str, layer: str) -> Dataset:
    dataset = DATASETS.get(f"{family}/{layer}")
    if dataset is None:
        fail(404, "unknown_dataset", f"No dataset {family}/{layer}")
    return dataset


DOWNLOAD_NAME = re.compile(r"^(?P<iso>[A-Z]{3})_adm(?P<level>[12])_(?P<slug>[a-z_]+)_metrics\.csv$")
DOWNLOAD_SLUGS = {d.replace("-", "_").replace("/", "_"): d for d in DATASETS}


def create_router(pool) -> tuple[APIRouter, Catalog]:
    """The upstream-compatible routes over a connection pool, and the catalog they keep in memory
    (call `catalog.refresh()` once the pool is open)."""
    router = APIRouter()
    catalog = Catalog(pool)

    def country(iso3: str) -> str:
        catalog.ensure_fresh()
        code = iso3.upper()
        if not any(c["iso3"] == code for c in catalog.countries):
            fail(404, "unknown_country", f"No data for {iso3}")
        return code

    def area(iso3: str, adm_id: str) -> tuple[str, int, str]:
        found = catalog.areas.get((queries.CGAZ, adm_id))
        if found is None or found[0] != country(iso3):
            fail(404, "unknown_area", f"No area {adm_id} in {iso3}")
        return found

    def interpret(dataset: Dataset, request: Request) -> tuple[dict, dict, dict]:
        query = dict(request.query_params)
        try:
            ctx, assignment = dataset.apply(query)
        except MappingError as error:
            fail(422, "invalid_filters", str(error))
        return query, ctx, assignment

    def series_values(series_id: int, iso3: str, level: int) -> list[tuple]:
        with pool.connection() as conn:
            return queries.series_values(conn, series_id, iso3, level)

    @router.get("/downloads/{iso3}/{filename}")
    def download(iso3: str, filename: str):
        """Card CSV, generated from the cube: one row per area x series.

        Same file names and key columns as the published research drop
        ({ISO3}_adm{level}_{slug}_metrics.csv), but dimension columns use the
        canonical names, not the published column layout.
        """
        code = country(iso3)
        match = DOWNLOAD_NAME.fullmatch(filename)
        dataset_id = match and match.group("iso") == code and DOWNLOAD_SLUGS.get(match.group("slug"))
        if not dataset_id:
            fail(404, "unknown_download", "unknown metrics download")
        level = int(match.group("level"))
        with pool.connection() as conn:
            rows = queries.download_rows(conn, dataset_id, code, level)
        dims = sorted({k for *_, d, _, _ in rows for k in d})
        name = next(c["name"] for c in catalog.countries if c["iso3"] == code)
        section, card = dataset_id.split("/")
        out = io.StringIO()
        writer = csv.writer(out)
        writer.writerow(["country_iso3", "country_name", "admin_level", "adm_id", "adm_name", "section", "card",
                         *dims, "metric", "unit", "value"])
        for adm_id, adm_name, metric, d, unit, value in rows:
            writer.writerow([code, name, f"ADM{level}", adm_id, adm_name, section, card,
                             *(d.get(k, "") for k in dims), metric, unit, "" if value is None else round(value, 3)])
        return Response(out.getvalue(), media_type="text/csv",
                        headers={"Content-Disposition": f'attachment; filename="{filename}"'})


    @router.get("/boundaries/{iso3}/{level}")
    def boundaries(iso3: str, level: int):
        code = country(iso3)
        key = (code, level)
        if key not in catalog.boundaries:
            with pool.connection() as conn:
                body = queries.boundaries_geojson(conn, code, level)
            if body is None:
                fail(404, "unknown_level", f"No admin level {level} for {iso3}")
            catalog.boundaries[key] = body
        return Response(catalog.boundaries[key], media_type="application/json",
                        headers={"Cache-Control": "public, max-age=3600"})


    @router.get(f"/metrics/{CURVES.path}/{{iso3}}")
    def concentration_curve(request: Request, iso3: str):
        code = country(iso3)
        query, ctx, assignment = interpret(CURVES, request)
        metric, dims = CURVES.canonical(assignment)
        series_id = catalog.series_id(CURVES.id, metric, dims)
        if series_id is None:
            fail(422, "invalid_filters", "No concentration curve for these parameters")
        with pool.connection() as conn:
            points = queries.curve_points(conn, series_id, code)
        return {**{k: query[k] for k in CURVE_ECHO if k in query},
                "points": [{"cumulative_population_share": x, "value": y} for x, y in points]}


    @router.get("/metrics/{family}/{layer}/{iso3}/profile/{adm_id}")
    def profile(request: Request, family: str, layer: str, iso3: str, adm_id: str):
        return breakdown(request, dataset_for(family, layer), iso3, adm_id, "profile")


    @router.get("/metrics/{family}/{layer}/{iso3}/ladder/{adm_id}")
    def ladder(request: Request, family: str, layer: str, iso3: str, adm_id: str):
        return breakdown(request, dataset_for(family, layer), iso3, adm_id, "ladder")


    @router.get("/metrics/{family}/{layer}/{iso3}/{level}")
    def series(request: Request, family: str, layer: str, iso3: str, level: int):
        dataset = dataset_for(family, layer)
        code = country(iso3)
        query, ctx, assignment = interpret(dataset, request)
        canonical = dataset.canonical({**assignment, **fill(dataset.series, ctx)})
        series_id = canonical and catalog.series_id(dataset.id, *canonical)
        if not series_id:
            fail(422, "invalid_filters", f"No {dataset.id} series for {query}")
        metric, dims = canonical
        rows = series_values(series_id, code, level)
        unit = catalog.units[(dataset.id, metric)]
        if ctx.get("$percent"):
            totals: dict[str, float] = {}
            for part in dataset.percent_of:
                part_id = catalog.series_id(dataset.id, metric, {**dims, **part})
                for adm, _, value in series_values(part_id, code, level) if part_id else []:
                    totals[adm] = totals.get(adm, 0.0) + (value or 0.0)
            rows = [(a, n, round(100 * v / totals[a], 3) if v is not None and totals.get(a) else None) for a, n, v in rows]
            unit = "percent"
        echo = {k: query[k] for k in dataset.echo.get("series", ()) if k in query}
        return {"unit": unit, **echo, "regions": [{"adm_id": a, "name": n, "value": v} for a, n, v in rows]}


    def breakdown(request: Request, dataset: Dataset, iso3: str, adm_id: str, kind: str):
        """Render one area's observations into the upstream profile or ladder layout."""
        _, _, name = area(iso3, adm_id)
        query, ctx, assignment = interpret(dataset, request)
        templates = dataset.profile if kind == "profile" else dataset.ladder
        if not templates:
            fail(404, "no_route", f"{dataset.id} has no {kind}")
        with pool.connection() as conn:
            rows = queries.area_observations(conn, dataset.id, adm_id)
        metric_order = {m: i for i, m in enumerate(dataset.metrics)}
        rows.sort(key=lambda r: (metric_order.get(r[0], 99), [catalog.order.get(kv, 999) for kv in sorted(r[1].items())]))
        targets = [{**{k: v for k, v in assignment.items() if k not in t.free(ctx)}, **t.fixed_for(ctx)} for t in templates]
        out: dict = {"adm_id": adm_id, "name": name}
        if kind == "ladder":
            out["unit"] = None
        out.update({k: query[k] for k in dataset.echo.get(kind, ()) if k in query})
        for template in templates:                           # literal top-level keys in template order
            first = template.path.split(".")[0]
            if "{" not in first:
                out.setdefault(first, None)
        for key in dataset.nulls:
            out.setdefault(key, None)
        for metric, dims, value in rows:
            for template, target in zip(templates, targets):
                if not dataset.matches(metric, dims, target):
                    continue
                keys = template.render({"metric": metric, **dims}, ctx)
                if keys is None:
                    continue
                unit = catalog.units[(dataset.id, metric)]
                node = out
                for key in keys[:-1]:
                    if not isinstance(node.get(key), dict):
                        node[key] = {}
                    node = node[key]
                if keys[-1] in node and node[keys[-1]] is not None:
                    continue
                if kind == "ladder":
                    node[keys[-1]] = value
                    out["unit"] = out["unit"] or unit
                else:
                    node[keys[-1]] = {"value": value, "unit": unit}
        if not any(key not in {"adm_id", "name", "unit", *query} and value is not None for key, value in out.items()):
            fail(422, "invalid_filters", f"No {dataset.id} {kind} values for {query}")
        return {k: v for k, v in out.items() if v is not None or k in dataset.nulls}

    return router, catalog
