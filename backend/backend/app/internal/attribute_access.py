from dataclasses import dataclass
from typing import Callable
from sqlalchemy import Column
from sqlalchemy.orm import Query
from sqlalchemy.sql import functions
from pydantic import Json, ValidationError


from backend.app import schemas
from backend.db import models


def add_damages_expected_value_query(
    fq: Query, dimensions: schemas.ExpectedDamagesDimensions, field: str
):
    q = fq.join(models.Feature.damages_expected)
    q = q.filter_by(rcp=dimensions.rcp, epoch=dimensions.epoch)
    agg = False
    if dimensions.hazard != "all":
        q = q.filter_by(hazard=dimensions.hazard)
    else:
        agg = True

    value: Column | functions.sum = getattr(models.ExpectedDamage, field)

    if agg:
        q = q.group_by(models.Feature.id)
        value = functions.sum(value)

    return q.add_column(value.label("value"))


# def add_damages_rp_value_query(fq: Query, dimesions: schemas.ReturnPeriodDamagesDimensions, field: str):
#     pass

# def add_damages_npv_value_query(fq: Query, dimesions: schemas.NPVDamagesDimensions, field: str):
#     pass

# def add_adaptation_value_query(fq: Query, dimesions: schemas.AdaptationDimensions, field: str):
#     pass


@dataclass
class DataGroupConfig:
    dimensions_schema: schemas.DataDimensions
    variables_schema: schemas.DataVariables
    add_value_query: Callable[[Query, schemas.DataDimensions, str], Query]


DATA_GROUP_CONFIGS: dict[str, DataGroupConfig] = {
    "damages_expected": DataGroupConfig(
        dimensions_schema=schemas.ExpectedDamagesDimensions,
        variables_schema=schemas.ExpectedDamagesVariables,
        add_value_query=add_damages_expected_value_query,
    ),
    # "damages_rp": DataGroupConfig(
    #     dimensions_schema=schemas.ReturnPeriodDamagesDimensions,
    #     variables_schema=schemas.ReturnPeriodDamagesVariables,
    #     add_value_query=add_damages_rp_value_query
    # ),
    # "damages_npv": DataGroupConfig(
    #     dimensions_schema=schemas.NPVDamagesDimensions,
    #     variables_schema=schemas.NPVDamagesVariables,
    #     add_value_query=add_damages_npv_value_query
    # ),
    # "adaptation": DataGroupConfig(
    #     dimensions_schema=schemas.AdaptationDimensions,
    #     variables_schema=schemas.AdaptationVariables,
    #     add_value_query=add_adaptation_value_query
    # )
}


def parse_dimensions(field_group: str, dimensions: Json):
    data_group_config = DATA_GROUP_CONFIGS.get(field_group)

    if data_group_config is not None:
        return data_group_config.dimensions_schema.parse_obj(dimensions)
    else:
        raise ValidationError(f"Invalid field group: {field_group}")


def add_value_query(
    q: Query, field_group: str, field_dimensions: schemas.DataDimensions, field: str
):
    data_group_config = DATA_GROUP_CONFIGS.get(field_group)

    if data_group_config is not None:
        return data_group_config.add_value_query(q, field_dimensions, field)
    else:
        raise ValidationError(f"Invalid field group: {field_group}")
