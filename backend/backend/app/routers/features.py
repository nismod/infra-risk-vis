from fastapi import APIRouter, Depends
from fastapi_pagination import Page, Params
from fastapi_pagination.ext.sqlalchemy import paginate
from sqlalchemy import desc
from sqlalchemy.orm import Session
from geoalchemy2 import functions

from backend.app import schemas
from backend.app.dependencies import get_db
from backend.app.internal.attribute_access import (
    add_value_query,
    parse_dimensions,
)
from backend.db import models


router = APIRouter(tags=["features"])


@router.get("/{feature_id}", response_model=schemas.FeatureOut)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    return feature


@router.get(
    "/sorted-by/{field_group}", response_model=Page[schemas.FeatureListItemOut[float]]
)
def read_sorted_features(
    layer: str,
    field_group: str,
    field: str,
    field_dimensions: schemas.DataDimensions = Depends(parse_dimensions),
    page_params: Params = Depends(),
    db: Session = Depends(get_db),
):
    base_query = (
        db.query(
            models.Feature.id.label("id"),
            models.Feature.string_id.label("string_id"),
            functions.ST_AsText(functions.Box2D(models.Feature.geom)).label("bbox_wkt"),
        )
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer)
    )

    q = add_value_query(base_query, field_group, field_dimensions, field).order_by(
        desc("value")
    )

    return paginate(q, page_params)
