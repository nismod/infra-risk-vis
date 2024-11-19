from dataclasses import dataclass
from typing import Callable, Optional
from sqlalchemy import Column, text, Float
from sqlalchemy.orm import Query
from sqlalchemy.sql import functions
from sqlalchemy.sql.operators import ColumnOperators
from pydantic import Json, ValidationError


from app import schemas
from db import models


def add_damages_expected_value_query(
    fq: Query,
    dimensions: schemas.ExpectedDamagesDimensions,
    field: str,
    field_params: schemas.DataParameters,
) -> Query:
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

import yaml

ADAPTATIONS_CONFIG = yaml.safe_load(
"""
properties:
  avoided_ead_amin:
    type: float
    json_key: "avoided_ead_amin"
  avoided_ead_mean:
    type: float
    json_key: "avoided_ead_mean"
  avoided_ead_amax:
    type: float
    json_key: "avoided_ead_amax"
  adaptation_cost:
    type: float
    json_key: "adaptation_cost"
  cost_benefit_ratio:
    type: calculated
    expression: "{avoided_ead_mean} / {adaptation_cost}"
"""
)


def build_sql_expression(data_column, data_config, field, params=None):
    """
    Build a SQL expression for the given field based on the configuration.
    """
    print(data_config)
    properties = data_config["properties"]
    field_config = properties.get(field)
    if not field_config:
        raise KeyError(f"Field '{field}' is not defined in configuration.")

    field_type = field_config["type"]
    if field_type == "calculated":
        # Substitute keys in the expression with their JSONB paths
        expression = field_config["expression"]
        for key, prop in properties.items():
            expression = expression.replace(
                f"{{{key}}}",
                f"({data_column[prop['json_key']]}.astext::float)"
            )
        # Substitute parameters if provided
        if params:
            for param_key, param_value in params.items():
                expression = expression.replace(f"{{{param_key}}}", str(param_value))
        return text(expression)
    elif field_type == "float":
        json_key = field_config["json_key"]
        return data_column[json_key].astext.cast(Float)
    else:
        raise ValueError(f"Unsupported field type '{field_type}'.")

def add_adaptation_value_query(
    fq: Query,
    dimensions: schemas.AdaptationDimensions,
    field: str,
    field_params: schemas.DataParameters,
):
    q = fq.join(models.Feature.adaptation)
    q = q.filter_by(
        hazard=dimensions.hazard,
        rcp=dimensions.rcp,
        adaptation_name=dimensions.adaptation_name,
        adaptation_protection_level=dimensions.adaptation_protection_level,
    )

    params = field_params.model_dump() if field_params else None
    value = build_sql_expression(models.AdaptationCostBenefit.properties, ADAPTATIONS_CONFIG, field, params)

    return q.add_column(value.label("value"))


@dataclass
class DataGroupConfig:
    dimensions_schema: schemas.DataDimensions
    variables_schema: schemas.DataVariables
    add_value_query: Callable[
        [Query, schemas.DataDimensions, str, schemas.DataParameters | None], Query
    ]
    field_parameters_schemas: Optional[dict[str, schemas.DataParameters]] = None


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
    "adaptation": DataGroupConfig(
        dimensions_schema=schemas.AdaptationDimensions,
        variables_schema=schemas.AdaptationVariables,
        add_value_query=add_adaptation_value_query,
        field_parameters_schemas={
            "cost_benefit_ratio": schemas.AdaptationCostBenefitRatioParameters,
        },
    ),
}


def parse_dimensions(field_group: str, dimensions: Json):
    data_group_config = DATA_GROUP_CONFIGS.get(field_group)

    if data_group_config is not None:
        return data_group_config.dimensions_schema.parse_obj(dimensions)
    else:
        raise ValidationError(f"Invalid field group: {field_group}")


def parse_parameters(field_group: str, field: str, parameters: Json):
    data_group_config = DATA_GROUP_CONFIGS.get(field_group)

    if data_group_config is not None:
        field_params_schema = data_group_config.field_parameters_schemas

        if field_params_schema is not None and field in field_params_schema:
            return field_params_schema[field].parse_obj(parameters)
        else:
            return None
    else:
        raise ValidationError(f"Invalid field group: {field_group}")


def add_value_query(
    q: Query,
    field_group: str,
    field_dimensions: schemas.DataDimensions,
    field: str,
    field_params: Optional[schemas.DataParameters] = None,
) -> Query:
    data_group_config = DATA_GROUP_CONFIGS.get(field_group)

    if data_group_config is not None:
        return data_group_config.add_value_query(
            q, field_dimensions, field, field_params
        )
    else:
        raise ValidationError(f"Invalid field group: {field_group}")
