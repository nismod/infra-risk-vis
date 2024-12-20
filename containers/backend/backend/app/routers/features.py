from typing import Optional

from fastapi import APIRouter, Depends, HTTPException
from fastapi_pagination import Page, Params
from fastapi_pagination.ext.sqlalchemy import paginate
from sqlalchemy import desc
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import Session
from geoalchemy2 import functions

from backend.app import schemas
from backend.app.dependencies import get_db
from backend.app.internal.attribute_access import (
    add_value_query,
    parse_dimensions,
    parse_parameters,
)
from backend.db import models


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
    layer: Optional[str] = None,
    sector: Optional[str] = None,
    subsector: Optional[str] = None,
    asset_type: Optional[str] = None,
):
    return schemas.LayerSpec(
        layer_name=layer,
        sector=sector,
        subsector=subsector,
        asset_type=asset_type,
    )


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
    db: Session = Depends(get_db),
):
    filled_layer_spec = {k: v for k, v in layer_spec.dict().items() if v is not None}
    base_query = (
        db.query(
            models.Feature.id.label("id"),
            models.Feature.string_id.label("string_id"),
            models.Feature.layer.label("layer"),
            functions.ST_AsText(functions.Box2D(models.Feature.geom)).label("bbox_wkt"),
        )
        .select_from(models.Feature)
        .join(models.FeatureLayer)
        .filter_by(**filled_layer_spec)
    )

    q = add_value_query(
        base_query, field_group, field_dimensions, field, field_params
    ).order_by(desc("value"))

    return paginate(q, page_params)
