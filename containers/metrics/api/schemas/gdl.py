"""
Pydantic Schemas for GDL metrics and boundaries
"""

from pydantic import BaseModel
from .geometry import MultiPolygon, Polygon

from typing import Literal

datasetKey = Literal["development", "education", "income", "healthcare"]


class AnnualData(BaseModel):
    """GDL dataset value for a region and year"""

    gdl_code: str
    region_name: str
    year: int
    value: float


class ValueExtent(BaseModel):
    """Min and and values for across dataset"""

    min: float
    max: float
    dataset: str


class CountryMeta(BaseModel):
    """GDL metadata for a country"""

    iso_code: str
    country_name: str
    continent: str


class RegionMeta(BaseModel):
    """GDL metadata for a region"""

    gdl_code: str
    region_name: str
    level: str
    iso_code: str


class NationalGeoProperties(RegionMeta):
    """Metadata to include with national geojson"""

    country_name: str


class NationalBoundary(BaseModel):
    """National geojson boundary feature"""

    type: str
    geometry: MultiPolygon


class NationalGeo(BaseModel):
    """National level geojson"""

    boundary: NationalBoundary
    envelope: Polygon
    properties: NationalGeoProperties


class SubnationalProperties(RegionMeta):
    """Metadata to include with sub-national geojson"""

    pass


class SubnationalGeo(BaseModel):
    """Region level geojson"""

    type: str
    geometry: MultiPolygon
    properties: SubnationalProperties
