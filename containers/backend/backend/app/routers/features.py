from typing import Any
from fastapi import APIRouter, Depends, HTTPException
from fastapi_pagination import Page, Params
from fastapi_pagination.ext.sqlalchemy import paginate
from pydantic import Json
from sqlalchemy import desc, Column, Integer
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import Session
from geoalchemy2 import functions

from app import schemas
from app.dependencies import get_db
from app.internal.attribute_access import (
    add_value_query,
    parse_dimensions,
    parse_parameters,
)
from db import models


router = APIRouter(tags=["features"])


@router.get("/{feature_id}", response_model=schemas.FeatureOut)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    try:
        feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    except NoResultFound:
        raise HTTPException(
            status_code=404,
            detail=f"No feature found with id {feature_id}",
        )
    return feature


def get_layer_spec(
    layer: str = None, sector: str = None, subsector: str = None, asset_type: str = None
):
    return schemas.LayerSpec(
        layer_name=layer,
        sector=sector,
        subsector=subsector,
        asset_type=asset_type,
    )

# parse json dictionary from ranking_scope parameter
def parse_ranking_scope(ranking_scope: Json = None):
    return ranking_scope or {}

def add_jsonb_filters(jsonb_column: Column, filters: dict[str, Any]):
    return [
        jsonb_column.op("->>")(key).cast(Integer) == value
        for key, value in filters.items()
    ]


@router.get(
    "/sorted-by/{field_group}", response_model=Page[schemas.FeatureListItemOut[float]]
)
def read_sorted_features(
    field_group: str,
    field: str,
    field_dimensions: schemas.DataDimensions = Depends(parse_dimensions),
    field_params: schemas.DataParameters = Depends(parse_parameters),
    layer_spec: schemas.LayerSpec = Depends(get_layer_spec),
    page_params: Params = Depends(),
    ranking_scope: dict = Depends(parse_ranking_scope),
    db: Session = Depends(get_db),
):
    filled_layer_spec = {k: v for k, v in layer_spec.model_dump().items() if v is not None}
    jsonb_filters = add_jsonb_filters(models.Feature.properties, ranking_scope)
    base_query = (
        db.query(
            models.Feature.id.label("id"),
            models.Feature.string_id.label("string_id"),
            models.Feature.layer.label("layer"),
            functions.ST_AsText(functions.Box2D(models.Feature.geom)).label("bbox_wkt"),
        )
        .select_from(models.Feature)
        .filter(*jsonb_filters)
        .join(models.FeatureLayer)
        .filter_by(**filled_layer_spec)
    )

    q = add_value_query(
        base_query, field_group, field_dimensions, field, field_params
    ).order_by(desc("value"))

    return paginate(q, page_params)
