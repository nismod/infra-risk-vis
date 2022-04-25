from typing import Any
from fastapi import APIRouter, Body, Depends
from sqlalchemy.orm import Session

from backend.app.dependencies import get_db
from backend.app.internal.attribute_access import (
    add_value_query,
    parse_dimensions,
)

from backend.db import models

from .. import schemas

router = APIRouter(tags=["attributes"])


@router.post("/{field_group}", response_model=schemas.AttributeLookup[Any | None])
def read_attributes(
    layer: str,
    field_group: str,
    field: str,
    field_dimensions: schemas.DataDimensions = Depends(parse_dimensions),
    ids: list[int] = Body(...),
    db: Session = Depends(get_db),
):
    base_query = (
        db.query(models.Feature.id)
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer, models.Feature.id.in_(ids))
    )
    query = add_value_query(base_query, field_group, field_dimensions, field)

    lookup = dict(query.all())

    return {id: lookup.get(id, None) for id in ids}
