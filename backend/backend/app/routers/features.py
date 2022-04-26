import json

from fastapi import APIRouter, Depends
from fastapi_pagination import Page, Params
from fastapi_pagination.ext.sqlalchemy import paginate
from pydantic import Json, ValidationError
from sqlalchemy import desc
from sqlalchemy.orm import Session, Query
from geoalchemy2 import functions

from backend.app import schemas
from backend.app.dependencies import get_damage_params, get_db
from backend.db import models


router = APIRouter(tags=["features"])


@router.get("/{feature_id}", response_model=schemas.FeatureOut)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    return feature


def get_field_params(field: str, field_params: Json):
    if field == "damages":
        print(field_params)
        return schemas.DamageParams.parse_obj(field_params)
    else:
        raise ValidationError(f"Invalid field: {field}")


def damages_query(q: Query, damage_params: schemas.DamageParams):
    return q.join(models.Feature.damages).filter_by(**damage_params.dict())


def get_field_query_params(field: str):
    if field == "damages":
        return models.Damage.mean, damages_query
    else:
        raise ValidationError(f"Invalid field: {field}")


@router.get(
    "/sorted-by/{field}", response_model=Page[schemas.FeatureListItemOut[float]]
)
def read_sorted_features(
    field: str,
    layer: str,
    field_params: schemas.FieldParams = Depends(get_field_params),
    page_params: Params = Depends(),
    db: Session = Depends(get_db),
):
    print(field_params) 
    value, augment_query = get_field_query_params(field)

    base_query = (
        db.query(
            models.Feature.id.label('id'),
            models.Feature.string_id.label('string_id'),
            functions.ST_AsText(functions.Box2D(models.Feature.geom)).label("bbox_wkt"),
            value.label("value"),
        )
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer)
    )
    print(base_query)
    q = augment_query(base_query, field_params).order_by(desc(value))

    print(q)
    return paginate(q, page_params)
