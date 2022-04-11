from fastapi import APIRouter, Body, Depends, Query
from sqlalchemy.orm import Session

from backend.app.dependencies import get_db
from backend.db import models

from .. import schemas

router = APIRouter()


def get_damage_params(
    hazard: str,
    rcp: str,
    epoch: str,
    damage_type: schemas.DamageType,
    protection_standard: int,
):
    args = locals()  # https://stackoverflow.com/a/2521937/1478817
    return schemas.DamageParams(**args)


@router.post("/damages", response_model=schemas.AttributeLookup[float | None])
def read_damages(
    layer: str,
    ids: list[int] = Body(...),
    params: schemas.DamageParams = Depends(get_damage_params),
    db: Session = Depends(get_db),
):
    lookup = dict(
        db.query(models.Feature.id, models.Damage.mean)
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer, models.Feature.id.in_(ids))
        .join(models.Feature.damages)
        .filter_by(**params.dict())
        .all()
    )

    return {id: lookup.get(id, None) for id in ids}
