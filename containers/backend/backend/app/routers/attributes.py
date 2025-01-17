from typing import Any
from fastapi import APIRouter, Body, Depends
from sqlalchemy import select

from backend.app.internal.attribute_access import (
    add_value_query,
    parse_dimensions,
    parse_parameters,
)

from backend.db import models
from backend.db.database import SessionDep

from backend.app import schemas

router = APIRouter(tags=["attributes"])


@router.post("/{field_group}", response_model=schemas.AttributeLookup[Any | None])
def read_attributes(
    layer: str,
    field_group: str,
    field: str,
    session: SessionDep,
    field_dimensions: schemas.DataDimensions = Depends(parse_dimensions),
    field_params: schemas.DataParameters = Depends(parse_parameters),
    ids: list[int] = Body(...),
):
    base_query = (
        select(models.Feature.id)
        .select_from(models.Feature)
        .filter(models.Feature.layer == layer, models.Feature.id.in_(ids))
    )
    query = add_value_query(
        base_query, field_group, field_dimensions, field, field_params
    )
    results = session.execute(query).all()

    lookup = dict(results)

    return {id: lookup.get(id, None) for id in ids}
