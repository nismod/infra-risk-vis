"""The /v1 API (contract C1) without HTTP: catalogue, coverage, and the four data reads.

The catalogue (`getCatalogue`) comes from the YAML dataset config
(`irv_metrics.config.load_catalogue`). The database supplies which countries
have data, the coverage of each, the series ids and the values. The cube is
seeded from the same config (`loaders.common.seed_catalogue`), so its dataset,
metric, dimension and member ids are the config's.

    service = Service(pool)            # anything with .connection() -> psycopg connection (search_path metrics)
    service.refresh()                  # after the pool is open; then ensure_fresh() before each request

Errors are raised as `ApiError(status, code, message)`; the router turns them
into `{"detail": {"code", "message"}}`.
"""
from __future__ import annotations

from collections import defaultdict
import csv
from dataclasses import dataclass, field
import hashlib
import io
import json
import logging
import threading
import time
from typing import Any

from .. import queries
from ..config import Catalogue, Dataset, load_catalogue

log = logging.getLogger(__name__)

RESERVED = frozenset({"level", "metric", "set", "area_id"})


class ApiError(Exception):
    def __init__(self, status: int, code: str, message: str):
        super().__init__(message)
        self.status, self.code, self.message = status, code, message


def not_found(code: str, message: str) -> ApiError:
    return ApiError(404, code, message)


def invalid(message: str) -> ApiError:
    return ApiError(422, "invalid_filters", message)


@dataclass(frozen=True)
class Series:
    id: int
    metric: str
    dims: dict[str, str]


@dataclass
class State:
    """Everything derived from one load version. Replaced whole on refresh, so a request sees one version."""

    version: str
    load_version: tuple
    catalogue_body: dict
    countries: dict[str, dict]                                   # iso3 -> CountryDetail body, countries with data
    set_kinds: dict[str, str]                                    # boundary set -> kind (config, then database)
    levels: dict[tuple[str, str], dict[int, int]]                # (iso3, set) -> {level: areas}
    series: dict[tuple[str, str, frozenset], int]                # (dataset, metric, dims) -> series_id
    by_dataset: dict[str, list[Series]]                          # dataset -> series, in catalogue order
    boundaries: dict[tuple, bytes] = field(default_factory=dict)  # (iso3, set, level) -> GeoJSON


class Service:
    """The /v1 reads over a connection pool.

    `countries`: restrict to these ISO3 codes (a single-country deployment); others are 404.
    `downloads`: whether `downloadDataset` answers (404 `downloads_disabled` when not).
    """

    def __init__(self, pool, catalogue: Catalogue | None = None, *, countries: set[str] | None = None,
                 downloads: bool = True, check_every: float = 5.0):
        self.pool = pool
        self.catalogue = catalogue or load_catalogue()
        self.only = {c.upper() for c in countries} if countries else None
        self.downloads = downloads
        self.check_every = check_every
        self.state: State | None = None
        self.checked = 0.0
        self._lock = threading.Lock()

    # ------------------------------------------------------------------ refresh

    def refresh(self):
        with self.pool.connection() as conn:
            load_version = tuple(queries.load_version(conn))
            state = self._build(conn, load_version)
        self.state = state
        self.checked = time.monotonic()

    def ensure_fresh(self):
        """Rebuild the state when a load has finished since the last check (checked every `check_every` s)."""
        if self.state is not None and time.monotonic() - self.checked < self.check_every:
            return
        with self._lock:
            if self.state is not None and time.monotonic() - self.checked < self.check_every:
                return
            with self.pool.connection() as conn:
                load_version = tuple(queries.load_version(conn))
                if self.state is None or load_version != self.state.load_version:
                    self.state = self._build(conn, load_version)
            self.checked = time.monotonic()

    def current(self) -> State:
        if self.state is None:
            self.ensure_fresh()
        return self.state

    def default_set(self) -> str | None:
        """The admin set for datasets whose source names none: the first admin set in the config."""
        return next((s.id for s in self.catalogue.boundary_sets.values() if s.kind == "admin"), None)

    def dataset_sets(self, dataset: Dataset) -> list[str]:
        default = self.default_set()
        out = []
        for s in dataset.boundary_sets(default) or [default]:
            if s is not None and s not in out:
                out.append(s)
        return out

    def _build(self, conn, load_version: tuple) -> State:
        cat = self.catalogue
        labels = queries.boundary_set_level_labels(conn)
        set_kinds = {s: kind for s, kind, _ in queries.boundary_set_rows(conn)}
        set_kinds.update({s.id: s.kind for s in cat.boundary_sets.values()})

        def level_label(boundary_set: str, level: int) -> str | None:
            spec = cat.boundary_sets.get(boundary_set)
            if spec is not None:
                for number, entry in zip(spec.level_numbers, spec.levels):
                    if number == level and entry.get("label"):
                        return entry["label"]
            return labels.get((boundary_set, level))

        body = self._catalogue_body(level_label)
        version = hashlib.sha256(json.dumps([body, [str(v) for v in load_version]], sort_keys=True,
                                            default=str).encode()).hexdigest()[:16]
        body = {"version": version, **body}

        # Series of the config's datasets and stored metrics (a cube seeded from another config may have more).
        series: dict[tuple, int] = {}
        by_dataset: dict[str, list[Series]] = defaultdict(list)
        for series_id, dataset_id, metric, dims in queries.series_list(conn):
            dataset = cat.datasets.get(dataset_id)
            if dataset is None or metric not in dataset.stored_metrics:
                continue
            series[(dataset_id, metric, frozenset(dims.items()))] = series_id
            by_dataset[dataset_id].append(Series(series_id, metric, dims))
        for dataset_id, items in by_dataset.items():
            items.sort(key=lambda s, d=cat.datasets[dataset_id]: self._series_order(d, s.metric, s.dims))

        levels: dict[tuple[str, str], dict[int, int]] = defaultdict(dict)
        for iso3, boundary_set, level, n in queries.country_set_levels(conn):
            levels[(iso3, boundary_set)][level] = n

        covered: dict[str, dict[str, dict[str, set[int]]]] = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
        for iso3, dataset_id, boundary_set, level in queries.coverage(conn):
            covered[iso3][dataset_id][boundary_set].add(level)

        countries = {}
        for iso3, name, admin_set, *bbox in queries.country_rows(conn):
            if self.only is not None and iso3 not in self.only:
                continue
            datasets = []
            for dataset in cat.datasets.values():
                found = covered.get(iso3, {}).get(dataset.id)
                if not found:
                    continue
                sets = [{"id": s, "levels": sorted(found[s])} for s in self.dataset_sets(dataset) if found.get(s)]
                if sets:
                    datasets.append({"id": dataset.id, "boundary_sets": sets, "download": self.downloads})
            if not datasets:
                continue                                        # countries with data only
            related = sorted((s for (i, s) in levels if i == iso3),
                             key=lambda s: (s != admin_set, list(cat.boundary_sets).index(s)
                                            if s in cat.boundary_sets else len(cat.boundary_sets), s))

            def summaries(boundary_set):
                return [{"level": lv, "label": level_label(boundary_set, lv), "areas": n}
                        for lv, n in sorted(levels[(iso3, boundary_set)].items())]

            countries[iso3] = {
                "iso3": iso3, "name": name,
                "bbox": tuple(bbox) if None not in bbox else (-180.0, -90.0, 180.0, 90.0),
                "boundary_set": admin_set,
                "levels": summaries(admin_set) if admin_set else [],
                "boundary_sets": [{"id": s, "levels": summaries(s)} for s in related],
                "datasets": datasets,
            }
        return State(version, load_version, body, countries, set_kinds, dict(levels), series, dict(by_dataset))

    def _catalogue_body(self, level_label) -> dict:
        cat = self.catalogue

        def compact_member(m):
            return {"id": m.id, "label": m.label, "short_label": m.short_label, "description": m.description,
                    "value": m.value, "code": m.code}

        datasets = []
        for d in cat.datasets.values():
            sets = self.dataset_sets(d)
            metrics = [{
                "id": m.id, "label": m.label, "unit": m.unit, "dimensions": list(d.metric_dimensions(m.id)),
                "aggregation": m.aggregation, "derived": m.derived, "aggregate_of": None,
                "boundary_sets": sets, "source_note": m.source_note, "description": m.description,
            } for m in d.metrics.values() if not (m.derive and "aggregate_of" in m.derive)]
            datasets.append({
                "id": d.id, "section": d.section, "title": d.title, "description": d.description,
                "documentation": d.documentation, "source": d.source_note, "kind": d.kind,
                "boundary_sets": sets, "metrics": metrics,
                "dimensions": [{"id": k, "members": list(v.members), "default": v.default, "labels": dict(v.labels)}
                               for k, v in d.dimensions.items()],
                "x": {"id": d.x["id"], "label": d.x["label"], "unit": d.x.get("unit")} if d.x else None,
                "display": dict(d.display) or None,
            })
        return {
            "sections": [{"id": s.id, "title": s.title, "description": s.description} for s in cat.sections.values()],
            "units": [{"id": u.id, "label": u.label, "quantity": u.quantity} for u in cat.units.values()],
            "dimensions": [{"id": dim.id, "label": dim.label, "description": dim.description, "ordered": dim.ordered,
                            "role": dim.role, "members": [compact_member(m) for m in dim.members.values()]}
                           for dim in cat.dimensions.values()],
            "boundary_sets": [{"id": s.id, "title": s.title, "kind": s.kind, "attribution": s.attribution,
                               "attribution_url": s.attribution_url,
                               "levels": [{"level": n, "label": level_label(s.id, n)} for n in s.level_numbers]}
                              for s in cat.boundary_sets.values()],
            "datasets": datasets,
        }

    @staticmethod
    def _series_order(dataset: Dataset, metric: str, dims: dict) -> tuple:
        """Metric in dataset order, then each of the metric's dimensions in member order (absent first)."""
        metrics = list(dataset.metrics)
        key = [metrics.index(metric) if metric in metrics else len(metrics)]
        for dim in dataset.metric_dimensions(metric):
            members = dataset.dimensions[dim].members if dim in dataset.dimensions else ()
            member = dims.get(dim)
            key.append(-1 if member is None else members.index(member) if member in members else len(members))
        return tuple(key)

    # ------------------------------------------------------------------ catalogue routes

    def list_countries(self) -> list[dict]:
        return [{k: c[k] for k in ("iso3", "name", "bbox", "boundary_set", "levels")}
                for c in self.current().countries.values()]

    def country(self, iso3: str) -> dict:
        found = self.current().countries.get(iso3.upper())
        if found is None:
            raise not_found("unknown_country", f"No data for country {iso3}")
        return found

    def catalogue_body(self) -> dict:
        return self.current().catalogue_body

    # ------------------------------------------------------------------ helpers

    def _dataset(self, section: str, name: str) -> Dataset:
        dataset = self.catalogue.datasets.get(f"{section}/{name}")
        if dataset is None:
            raise not_found("unknown_dataset", f"No dataset {section}/{name}")
        return dataset

    def _set(self, dataset: Dataset, boundary_set: str | None) -> str:
        sets = self.dataset_sets(dataset)
        if boundary_set is None:
            return sets[0]
        if boundary_set not in sets:
            raise not_found("unknown_boundary_set", f"{dataset.id} has no data on boundary set {boundary_set}"
                                                    f" (it has {', '.join(sets)})")
        return boundary_set

    def _level(self, state: State, iso3: str, boundary_set: str, level: int):
        if level not in state.levels.get((iso3, boundary_set), {}):
            raise not_found("unknown_level", f"No level {level} of {boundary_set} in {iso3}")

    def _share(self, state: State, boundary_set: str, share: float | None) -> float | None:
        return share if state.set_kinds.get(boundary_set) == "overlay" else None

    @staticmethod
    def _dims(dataset: Dataset, query: list[tuple[str, str]]) -> dict[str, str]:
        """Dimension filters from the query string: every parameter that is not reserved."""
        dims: dict[str, str] = {}
        for key, value in query:
            if key in RESERVED:
                continue
            if key not in dataset.dimensions:
                known = ", ".join(dataset.dimensions) or "none"
                raise invalid(f"{dataset.id} has no dimension {key} (dimensions: {known})")
            if key in dims:
                raise invalid(f"{key} is given more than once")
            if value not in dataset.dimensions[key].members:
                raise invalid(f"{key}={value} is not a member for {dataset.id}"
                              f" (members: {', '.join(dataset.dimensions[key].members)})")
            dims[key] = value
        return dims

    @staticmethod
    def _metric(dataset: Dataset, metric: str | None) -> str:
        if metric is None:
            stored = list(dataset.stored_metrics)
            if len(stored) != 1:
                raise invalid(f"{dataset.id} has several metrics; give metric= ({', '.join(stored)})")
            return stored[0]
        spec = dataset.metrics.get(metric)
        if spec is None or (spec.derive and "share_of" not in spec.derive):
            known = [m for m, s in dataset.metrics.items() if not (s.derive and "share_of" not in s.derive)]
            raise invalid(f"{dataset.id} has no metric {metric} (metrics: {', '.join(known)})")
        return metric

    @staticmethod
    def _keyed(dataset: Dataset, metric: str, dims: dict) -> dict:
        keyed = dataset.metric_dimensions(metric)
        return {k: v for k, v in dims.items() if k in keyed}

    def _series_id(self, state: State, dataset: Dataset, metric: str, dims: dict) -> int | None:
        return state.series.get((dataset.id, metric, frozenset(dims.items())))

    def _one_series(self, state: State, dataset: Dataset, metric: str, dims: dict) -> tuple[str, dict, int]:
        """(stored metric, the series' dims, series_id) for a complete assignment; 422 when it names no series."""
        spec = dataset.metrics[metric]
        base = spec.derive["share_of"] if spec.derive else metric
        dims = self._keyed(dataset, base, dims)
        series_id = self._series_id(state, dataset, base, dims)
        if series_id is None:
            keyed = dataset.metric_dimensions(base)
            missing = [d for d in keyed if d not in dims]
            hint = f"; the series of {metric} may be keyed by {', '.join(missing)}" if missing else ""
            raise invalid(f"No {dataset.id} series for metric={metric} {json.dumps(dims, sort_keys=True)}{hint}")
        return base, dims, series_id

    def _totals(self, state: State, dataset: Dataset, metric: str, dims: dict) -> list[int]:
        """The series whose values sum to the denominator of a share_of metric, for one base series."""
        derive = dataset.metrics[metric].derive
        out = []
        for part in derive["total"]:
            series_id = self._series_id(state, dataset, derive["share_of"], {**dims, **part})
            if series_id is not None:
                out.append(series_id)
        return out

    @staticmethod
    def _ratio(dataset: Dataset, metric: str, value, total) -> float | None:
        if value is None or not total:
            return None
        return (100.0 if dataset.metrics[metric].unit == "percent" else 1.0) * value / total

    def _scope(self, iso3: str, section: str, name: str, boundary_set: str | None, kind: str):
        state = self.current()
        code = self.country(iso3)["iso3"]
        dataset = self._dataset(section, name)
        if dataset.kind != kind:
            route = "curve" if dataset.kind == "curve" else "values, areas"
            raise not_found("unknown_dataset", f"{dataset.id} is a {dataset.kind} dataset; use {route}")
        return state, code, dataset, self._set(dataset, boundary_set)

    def _areas(self, state: State, code: str, boundary_set: str, area_ids: list[str]) -> dict[str, tuple]:
        with self.pool.connection() as conn:
            rows = queries.areas_in_country(conn, boundary_set, area_ids, code)
        found = {r[0]: r for r in rows}
        missing = [a for a in dict.fromkeys(area_ids) if a not in found]
        if missing:
            raise not_found("unknown_area", f"No area {', '.join(missing)} of {boundary_set} in {code}")
        return found

    # ------------------------------------------------------------------ data routes

    def values(self, iso3, section, name, level, boundary_set, metric, query) -> dict:
        state, code, dataset, boundary_set = self._scope(iso3, section, name, boundary_set, "values")
        dims = self._dims(dataset, query)
        metric = self._metric(dataset, metric)
        base, dims, series_id = self._one_series(state, dataset, metric, dims)
        self._level(state, code, boundary_set, level)
        derived = metric != base
        totals = self._totals(state, dataset, metric, dims) if derived else []
        with self.pool.connection() as conn:
            rows = queries.series_area_values(conn, [series_id, *totals], code, boundary_set, level)
        by_area: dict[str, dict] = {}
        for area_id, area_name, share, sid, value in rows:
            entry = by_area.setdefault(area_id, {"area_id": area_id, "name": area_name,
                                                 "share": self._share(state, boundary_set, share), "v": {}})
            entry["v"][sid] = value
        values = []
        for entry in by_area.values():
            found = entry.pop("v")
            value = found.get(series_id)
            if derived:
                value = self._ratio(dataset, metric, value, sum(found.get(t) or 0.0 for t in totals))
            values.append({**entry, "value": value})
        return {"iso3": code, "dataset": dataset.id, "metric": metric, "dims": dims,
                "unit": dataset.metrics[metric].unit, "level": level, "boundary_set": boundary_set, "values": values}

    def area_observations(self, iso3, section, name, area_ids, boundary_set, metric, query) -> dict:
        state, code, dataset, boundary_set = self._scope(iso3, section, name, boundary_set, "values")
        dims = self._dims(dataset, query)
        if metric is not None:
            metric = self._metric(dataset, metric)
        areas = self._areas(state, code, boundary_set, area_ids)
        spec = dataset.metrics.get(metric) if metric else None
        base = spec.derive["share_of"] if spec is not None and spec.derive else metric
        wanted = [s for s in state.by_dataset.get(dataset.id, ())
                  if (base is None or s.metric == base)
                  and all(s.dims.get(k) == v for k, v in self._keyed(dataset, s.metric, dims).items())]
        totals = {s.id: self._totals(state, dataset, metric, s.dims) for s in wanted} if base != metric else {}
        fetch = {s.id for s in wanted} | {t for ids in totals.values() for t in ids}
        stored: dict[tuple[int, int], Any] = {}
        if fetch:
            with self.pool.connection() as conn:
                for area_key, series_id, value in queries.area_series_values(
                        conn, [r[1] for r in areas.values()], sorted(fetch)):
                    stored[(area_key, series_id)] = value
        entries = []
        for area_id in area_ids:
            _, area_key, level, area_name, share = areas[area_id]
            observations = []
            for s in wanted:
                if (area_key, s.id) not in stored:
                    continue
                value = stored[(area_key, s.id)]
                if s.id in totals:
                    value = self._ratio(dataset, metric, value,
                                        sum(stored.get((area_key, t)) or 0.0 for t in totals[s.id]))
                shown = metric if s.id in totals else s.metric
                observations.append({"metric": shown, "dims": s.dims, "value": value,
                                     "unit": dataset.metrics[shown].unit})
            entries.append({"area_id": area_id, "name": area_name, "level": level,
                            "share": self._share(state, boundary_set, share), "observations": observations})
        return {"iso3": code, "dataset": dataset.id, "boundary_set": boundary_set, "areas": entries}

    def curve(self, iso3, section, name, boundary_set, metric, area_id, query) -> dict:
        state, code, dataset, boundary_set = self._scope(iso3, section, name, boundary_set, "curve")
        dims = self._dims(dataset, query)
        metric = self._metric(dataset, metric)
        if area_id is None:
            if state.set_kinds.get(boundary_set) != "admin":
                raise invalid(f"area_id is required on {boundary_set}, which is not an admin set")
            area_id = code
        _, dims, series_id = self._one_series(state, dataset, metric, dims)
        _, area_key, *_ = self._areas(state, code, boundary_set, [area_id])[area_id]
        with self.pool.connection() as conn:
            points = queries.curve(conn, series_id, area_key)
        x = dataset.x or {"id": "x", "label": "x"}
        return {"iso3": code, "dataset": dataset.id, "metric": metric, "dims": dims, "area_id": area_id,
                "boundary_set": boundary_set, "unit": dataset.metrics[metric].unit,
                "x": {"id": x["id"], "label": x["label"], "unit": x.get("unit")},
                "points": [{"x": px, "y": py} for px, py in points]}

    def boundaries(self, iso3: str, level: int, boundary_set: str | None) -> bytes:
        state = self.current()
        country = self.country(iso3)
        code = country["iso3"]
        boundary_set = boundary_set or country["boundary_set"]
        if not state.levels.get((code, boundary_set)):
            raise not_found("unknown_boundary_set", f"No areas of boundary set {boundary_set} in {code}")
        if state.set_kinds.get(boundary_set) == "group":
            raise not_found("unknown_boundary_set", f"{boundary_set} is a group set; it has no geometry")
        self._level(state, code, boundary_set, level)
        key = (code, boundary_set, level)
        body = state.boundaries.get(key)
        if body is None:
            with self.pool.connection() as conn:
                body = queries.area_features(conn, code, boundary_set, level,
                                             with_share=state.set_kinds.get(boundary_set) == "overlay")
            if body is None:
                raise not_found("unknown_level", f"No level {level} of {boundary_set} in {code}")
            state.boundaries[key] = body
        return body

    def download(self, iso3, section, name, level, boundary_set) -> tuple[str, str]:
        """(file name, CSV text) of a dataset for one country, set and level."""
        if not self.downloads:
            raise not_found("downloads_disabled", "Downloads are turned off in this deployment")
        state = self.current()
        country = self.country(iso3)
        code = country["iso3"]
        dataset = self._dataset(section, name)
        boundary_set = self._set(dataset, boundary_set)
        self._level(state, code, boundary_set, level)
        series = {s.id: s for s in state.by_dataset.get(dataset.id, ())}
        order = {sid: i for i, sid in enumerate(series)}
        overlay = state.set_kinds.get(boundary_set) == "overlay"
        dims = [d for d in dataset.dimensions if any(d in s.dims for s in series.values())]
        curve = dataset.kind == "curve"
        x_id = (dataset.x or {}).get("id", "x")
        with self.pool.connection() as conn:
            if curve:
                rows = queries.level_curve_points(conn, list(series), code, boundary_set, level)
            else:
                rows = [(*r, None) for r in queries.level_observations(conn, list(series), code, boundary_set, level)]
        if not curve:
            rows.sort(key=lambda r: (r[0], order[r[3]]))
        out = io.StringIO()
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["country_iso3", "boundary_set", "level", "area_id", "area_name",
                         *(["share"] if overlay else []), *dims, "metric", "unit",
                         *([x_id] if curve else []), "value"])
        for row in rows:
            area_id, area_name, share, series_id = row[:4]
            s = series[series_id]
            values = [row[4], row[5]] if curve else [row[4]]
            writer.writerow([code, boundary_set, level, area_id, area_name, *([share] if overlay else []),
                             *(s.dims.get(d, "") for d in dims), s.metric, dataset.metrics[s.metric].unit,
                             *("" if v is None else repr(v) for v in values)])
        filename = f"{code}_{boundary_set}_level{level}_{dataset.id.replace('/', '_')}.csv"
        return filename, out.getvalue()
