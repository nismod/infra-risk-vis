"""Response models and parameter types of the /v1 API (contract C1).

A copy of the declarations in `backend/contract/v1.py`, which is not part of
this package (GRI vendors `irv_metrics` without `contract/`). Keep them
identical: `tests/test_contract_conformance.py` builds the OpenAPI schema of
`irv_metrics.api.create_router` and compares it with `backend/openapi.json`.
"""
from __future__ import annotations

from typing import Annotated, Any, Literal

from fastapi import Path, Query
from pydantic import BaseModel, Field

Iso3 = Annotated[str, Path(pattern="^[A-Za-z]{3}$", description="ISO 3166-1 alpha-3 code, any case.")]
Section = Annotated[str, Path(description="Section id, e.g. `hazard`.")]
DatasetName = Annotated[str, Path(description="Dataset name within the section, e.g. `river-flooding`.")]
Level = Annotated[int, Query(ge=0, description="Level within the boundary set. In admin sets, 0 is the country.")]
BoundarySet = Annotated[str | None, Query(
    alias="set", description="Boundary set id. Defaults to the dataset's first set (data routes) or the "
                             "country's admin set (`getBoundaries`).")]
AreaIds = Annotated[list[str], Query(
    min_length=1, max_length=50,
    description="Area ids in the chosen set, repeated: `?area_id=KEN&area_id=KEN.1_1`. Each must be related "
                "to the country. The country itself is `area_id={ISO3}` in an admin set.")]
Metric = Annotated[str | None, Query(description="Metric id. May be omitted when the dataset stores one metric.")]

DIMS = {
    "name": "dims", "in": "query", "required": False, "style": "form", "explode": True,
    "description": "Dimension filters, one query parameter per dimension: `?sector=total&return_period_years=100`.",
    "schema": {"type": "object", "additionalProperties": {"type": "string"}},
}


# --------------------------------------------------------------------------- models

class ErrorDetail(BaseModel):
    code: Literal["unknown_country", "unknown_dataset", "unknown_area", "unknown_level", "unknown_boundary_set",
                  "invalid_filters", "downloads_disabled", "busy"]
    message: str


class Error(BaseModel):
    detail: ErrorDetail


class LevelSummary(BaseModel):
    level: int
    label: str | None = Field(None, description="Local name of the level, e.g. `County`, `Parish`.")
    areas: int = Field(description="Areas at this level related to the country.")


class Country(BaseModel):
    iso3: str = Field(description="Upper case.")
    name: str
    bbox: tuple[float, float, float, float] = Field(description="west, south, east, north.")
    boundary_set: str = Field(description="The country's admin set: the default for maps and boundaries.")
    levels: list[LevelSummary] = Field(description="Levels of the admin set, 0 first.")


class CountryBoundarySet(BaseModel):
    id: str
    levels: list[LevelSummary]


class SetCoverage(BaseModel):
    id: str = Field(description="Boundary set id.")
    levels: list[int] = Field(description="Levels with data.")


class DatasetCoverage(BaseModel):
    id: str = Field(description="`{section}/{dataset}`.")
    boundary_sets: list[SetCoverage] = Field(description="Sets with data for this country, in the dataset's order.")
    download: bool = Field(description="Whether the download route answers for this dataset.")


class CountryDetail(Country):
    boundary_sets: list[CountryBoundarySet] = Field(
        description="Every set with areas related to this country, the admin set first.")
    datasets: list[DatasetCoverage] = Field(description="Datasets with data for this country. Others are absent.")


class SectionModel(BaseModel):
    id: str
    title: str
    description: str | None = None


class Unit(BaseModel):
    id: str
    label: str = Field(description="Display form, e.g. `km²`.")
    quantity: str = Field(description="e.g. `area`, `count`, `currency`, `ratio`.")


class Member(BaseModel):
    id: str
    label: str
    short_label: str | None = None
    description: str | None = None
    value: float | None = Field(None, description="Numeric value of the member, e.g. 100 for a 100-year return period.")
    code: str | None = Field(None, description="External code, e.g. a GHS-SMOD class.")


class Dimension(BaseModel):
    id: str
    label: str
    description: str | None = None
    ordered: bool = Field(description="Members have a meaningful order (return periods, years, quintiles).")
    role: Literal["time"] | None = None
    members: list[Member] = Field(description="All members, in display order.")


class MetricModel(BaseModel):
    id: str
    label: str
    unit: str
    dimensions: list[str] = Field(description="The dimensions that key this metric.")
    aggregation: Literal["sum", "mean", "weighted_mean", "none"]
    derived: bool = Field(False, description="Computed on read. Returned only when asked for by `metric=`.")
    aggregate_of: str | None = Field(
        None, description="Set when this metric's values are aggregated at load from another metric of the dataset "
                          "over the member countries of each group, e.g. a population-weighted mean of national "
                          "values. Published group figures are separate metrics, without this field.")
    boundary_sets: list[str] = Field(description="Sets this metric has data on, in the dataset's order.")
    source_note: str | None = Field(None, description="Where this metric's values come from, when it differs "
                                                      "from the dataset's `source`.")
    description: str | None = None


class DatasetDimension(BaseModel):
    id: str
    members: list[str] = Field(description="The members this dataset uses, in display order.")
    default: str | None = None
    labels: dict[str, str] = Field(default_factory=dict, description="Label overrides for this dataset, by member id.")


class CurveAxis(BaseModel):
    id: str
    label: str
    unit: str | None = None


class Dataset(BaseModel):
    id: str = Field(description="`{section}/{dataset}`.")
    section: str
    title: str
    description: str | None = None
    documentation: str | None = Field(None, description="URL or host-relative path of the method notes.")
    source: str | None = None
    kind: Literal["values", "curve"] = Field(
        description="`values`: one number per series and area. `curve`: a list of points per series and area.")
    boundary_sets: list[str] = Field(
        description="Sets this dataset has data on: the union of its metrics' sets. The first is the default "
                    "for `set`.")
    metrics: list[MetricModel]
    dimensions: list[DatasetDimension]
    x: CurveAxis | None = Field(None, description="`curve` datasets only.")
    display: dict[str, Any] | None = Field(
        None, description="Display hints, as written in the dataset config (contract C2). Opaque to the API.")


class BoundarySetLevel(BaseModel):
    level: int
    label: str | None = None


class BoundarySetModel(BaseModel):
    id: str = Field(description="e.g. `cgaz-6.0.0`, `hydrobasins-v1c`.")
    title: str
    kind: Literal["admin", "overlay", "group"] = Field(description=(
        "`admin`: nested areas that partition each country; level 0 is the country, with `area_id` = ISO3. "
        "`overlay`: areas drawn independently of borders (river basins); one area may relate to several "
        "countries, and the areas related to a country do not partition it, so values must not be summed "
        "across them. `group`: named sets of countries (continents, income groups), without geometry of "
        "their own."))
    attribution: str | None = None
    attribution_url: str | None = None
    levels: list[BoundarySetLevel]


class Catalogue(BaseModel):
    version: str = Field(description="Changes whenever data or config is reloaded. Also the ETag.")
    sections: list[SectionModel]
    units: list[Unit]
    dimensions: list[Dimension]
    boundary_sets: list[BoundarySetModel]
    datasets: list[Dataset]


SHARE = Field(None, description="Overlay sets: the fraction of the area that lies in this country, 0–1. "
                                "Null in admin and group sets.")


class AreaValue(BaseModel):
    area_id: str
    name: str
    value: float | None
    share: float | None = SHARE


class SeriesValues(BaseModel):
    iso3: str = Field(description="The country the request is scoped to.")
    dataset: str
    metric: str
    dims: dict[str, str] = Field(description="The full dimension assignment of the series.")
    unit: str
    level: int
    boundary_set: str
    values: list[AreaValue] = Field(
        description="Every area at this level related to the country, ordered by `area_id`. Overlay areas "
                    "carry their whole-area value, not the part inside the country.")


class Observation(BaseModel):
    metric: str
    dims: dict[str, str]
    value: float | None
    unit: str


class AreaEntry(BaseModel):
    area_id: str
    name: str
    level: int
    share: float | None = SHARE
    observations: list[Observation] = Field(
        description="Ordered by metric (dataset order), then by each dimension's member order.")


class AreaObservations(BaseModel):
    iso3: str = Field(description="The country the request is scoped to.")
    dataset: str
    boundary_set: str
    areas: list[AreaEntry] = Field(description="One entry per requested `area_id`, in request order.")


class CurvePoint(BaseModel):
    x: float
    y: float


class Curve(BaseModel):
    iso3: str
    dataset: str
    metric: str
    dims: dict[str, str]
    area_id: str
    boundary_set: str
    unit: str
    x: CurveAxis
    points: list[CurvePoint]


NOT_FOUND = {404: {"model": Error}}
FILTERED = {404: {"model": Error}, 422: {"model": Error}}
