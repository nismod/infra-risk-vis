import json
from typing import List
from fastapi import APIRouter, HTTPException
import logging
from api.database.database import get_db
from api.database import gdl as models
from sqlalchemy.sql import select, func
from api.schemas import gdl as schemas
from sqlalchemy.orm import Session
from fastapi import APIRouter, Depends, HTTPException

# Retrieval of all country metrics data
API_ROUTE_BASE = "/metrics"

# GDL geojson
GDL_BOUNDARY_BASE_ROUTE = API_ROUTE_BASE + "/geojson"
GDL_SUBNATIONAL_FROM_ISO_ROUTE = GDL_BOUNDARY_BASE_ROUTE + "/subnational/iso/{iso_code}"
GDL_NATIONAL_FROM_ISO_ROUTE = GDL_BOUNDARY_BASE_ROUTE + "/national/iso/{iso_code}"

# GDL countries and regions metadata
GDL_META_BASE_ROUTE = API_ROUTE_BASE + "/meta"
GDL_COUNTRY_META_ROUTE = GDL_META_BASE_ROUTE + "/countries"
GDL_REGION_META_ROUTE = GDL_META_BASE_ROUTE + "/regions"

# GDL annual metrics
GDL_DATA_BASE_ROUTE = API_ROUTE_BASE + "/data"
GDL_DATA_ISO_ROUTE = GDL_DATA_BASE_ROUTE + "/{dataset_key}" + "/{iso_code}"
GDL_DATA_EXTENT_ROUTE = GDL_DATA_BASE_ROUTE + "/{dataset_key}" + "/extent"


router = APIRouter(
    tags=["boundaries"],
    dependencies=[],
    responses={404: {"description": "Not found"}},
)


def parse_gdl_annual(results, dataset_key):
    return [
        {
            "gdl_code": result.gdl_code,
            "region_name": result.region_name,
            "year": result.year,
            "value": getattr(result, dataset_key),
        }
        for result in results
    ]


def select_gdl_annual_join_region():
    return select(
        models.GdlAnnual.gdl_code,
        models.GdlAnnual.year,
        models.GdlAnnual.development,
        models.GdlAnnual.education,
        models.GdlAnnual.income,
        models.GdlAnnual.healthcare,
        models.GdlRegion.region_name,
    ).join(models.GdlRegion, models.GdlAnnual.gdl_code == models.GdlRegion.gdl_code)


@router.get(GDL_DATA_EXTENT_ROUTE, response_model=schemas.ValueExtent)
async def get_data_extent(
    dataset_key: schemas.datasetKey, db: Session = Depends(get_db)
):
    """Calculate value extent across all countries and years for a single dataset"""

    stmt = select_gdl_annual_join_region()
    results = db.execute(stmt)

    # pool non-null values for dataset
    values = [getattr(result, dataset_key) for result in results]
    filtered_values = list(filter((lambda x: isinstance(x, float)), values))

    minimum = min(filtered_values)
    maximum = max(filtered_values)

    return {"min": minimum, "max": maximum, "dataset": dataset_key}


@router.get(GDL_DATA_ISO_ROUTE, response_model=List[schemas.AnnualData])
async def get_annual_for_country(
    dataset_key: schemas.datasetKey, iso_code: str, db: Session = Depends(get_db)
):
    """Read all annual data for single country and dataset"""
    try:
        stmt = (
            select_gdl_annual_join_region()
            .where(models.GdlRegion.iso_code == iso_code)
            .order_by(models.GdlAnnual.gdl_code, models.GdlAnnual.year)
        )
        results = db.execute(stmt)

        data = [
            {
                "gdl_code": result.gdl_code,
                "region_name": result.region_name,
                "year": result.year,
                "value": getattr(result, dataset_key),
            }
            for result in results
        ]

        return data

    except Exception as err:
        raise HTTPException(status_code=500)


@router.get(GDL_COUNTRY_META_ROUTE, response_model=List[schemas.CountryMeta])
async def get_all_countries_meta(db: Session = Depends(get_db)):
    """Read all countries metadata"""
    try:
        stmt = select(models.IsoCountry)
        results = db.execute(stmt).scalars()

        data = [
            {
                "iso_code": result.iso_code,
                "country_name": result.country_name,
                "continent": result.continent,
            }
            for result in results
        ]

        return data

    except Exception as err:
        raise HTTPException(status_code=500)


@router.get(GDL_REGION_META_ROUTE, response_model=List[schemas.RegionMeta])
async def get_all_regions_meta(db: Session = Depends(get_db)):
    """Read all region metadata"""
    try:
        stmt = select(models.GdlRegion)
        results = db.execute(stmt).scalars()

        data = [
            {
                "gdl_code": result.gdl_code,
                "region_name": result.region_name,
                "level": result.level,
                "iso_code": result.iso_code,
            }
            for result in results
        ]

        return data

    except Exception as err:
        raise HTTPException(status_code=500)


@router.get(GDL_NATIONAL_FROM_ISO_ROUTE, response_model=schemas.NationalGeo)
async def get_national_for_iso(iso_code: str, db: Session = Depends(get_db)):
    """Read all national GeoJSON tied to a single country iso_code"""
    try:
        stmt = (
            select(
                func.ST_AsGeoJSON(models.GdlNational.geometry),
                func.ST_AsGeoJSON(func.ST_Envelope(models.GdlNational.geometry)),
                models.GdlNational.gdl_code,
                models.GdlRegion.region_name,
                models.GdlRegion.level,
                models.GdlRegion.iso_code,
                models.IsoCountry.country_name,
            )
            .join(
                models.GdlRegion,
                models.GdlRegion.gdl_code == models.GdlNational.gdl_code,
            )
            .join(
                models.IsoCountry,
                models.GdlRegion.iso_code == models.IsoCountry.iso_code,
            )
            .where(models.GdlRegion.iso_code == iso_code)
        )

        result = db.execute(stmt).first()

        if result is None:
            raise HTTPException(status_code=404, detail="Boundary not found")

        data = {
            "boundary": {
                "type": "Feature",
                "geometry": json.loads(result[0]),
            },
            "envelope": json.loads(result[1]),
            "properties": {
                "gdl_code": result.gdl_code,
                "region_name": result.region_name,
                "level": result.level,
                "iso_code": result.iso_code,
                "country_name": result.country_name,
            },
        }

        return data

    except HTTPException as not_found:
        raise not_found

    except Exception as err:
        raise HTTPException(status_code=500)


@router.get(GDL_SUBNATIONAL_FROM_ISO_ROUTE, response_model=List[schemas.SubnationalGeo])
async def get_all_boundaries_for_iso(iso_code: str, db: Session = Depends(get_db)):
    """Read all GeoJSON tied to a single country iso_code"""
    try:
        stmt = (
            select(
                func.ST_AsGeoJSON(models.GdlSubnational.geometry),
                models.GdlSubnational.gdl_code,
                models.GdlRegion.region_name,
                models.GdlRegion.level,
                models.GdlRegion.iso_code,
            )
            .join(models.GdlRegion)
            .where(models.GdlRegion.iso_code == iso_code)
        )
        results = db.execute(stmt)

        if results is None:
            raise HTTPException(status_code=404, detail="No boundaries found")

        data = [
            {
                "type": "Feature",
                "geometry": json.loads(result[0]),
                "properties": {
                    "gdl_code": result.gdl_code,
                    "region_name": result.region_name,
                    "level": result.level,
                    "iso_code": result.iso_code,
                },
            }
            for result in results
        ]

        return data

    except HTTPException as not_found:
        raise not_found

    except Exception as err:
        raise HTTPException(status_code=500)
