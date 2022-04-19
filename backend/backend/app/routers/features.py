from fastapi import APIRouter, Depends
from fastapi_pagination import Page, Params
from fastapi_pagination.ext.sqlalchemy import paginate
from sqlalchemy import desc
from sqlalchemy.orm import Session
from geoalchemy2 import functions

from backend.app import schemas
from backend.app.dependencies import get_damage_params, get_db
from backend.db import models


router = APIRouter(tags=["features"])


@router.get("/{feature_id}", response_model=schemas.FeatureOut)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    return feature


@router.get(
    "/sorted-by/damages", response_model=Page[schemas.FeatureListItemOut[float]]
)
def read_sorted_features(
    layer: str,
    damage_params: schemas.DamageParams = Depends(get_damage_params),
    page_params: Params = Depends(),
    db: Session = Depends(get_db),
):
    q = (
        db.query(
            models.Feature.id,
            models.Feature.string_id,
            functions.ST_AsText(functions.Box2D(models.Feature.geom)).label("bbox_wkt"),
            models.Damage.mean.label("value"),
        )
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer)
        .join(
            models.Feature.damages,
        )
        .filter_by(**damage_params.dict())
        .order_by(desc("value"))
    )

    return paginate(q, page_params)
