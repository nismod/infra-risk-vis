from typing import Any, Literal
from pydantic import BaseModel
from sqlalchemy import Column, Float, cast, literal_column

class StaticPropertyConfig(BaseModel):
    type: str
    json_key: str

class CalculatedPropertyConfig(BaseModel):
    type: Literal['calculated']
    expression: str

class DynamicDataConfig(BaseModel):
    properties: dict[str, StaticPropertyConfig | CalculatedPropertyConfig]

def build_sql_expression(data_column: Column, data_config: DynamicDataConfig, field, params: dict[str, Any]=None):
    """
    Build a SQL expression for the given field based on the configuration.
    """
    properties_config = data_config.properties
    field_config = properties_config.get(field)
    if not field_config:
        raise KeyError(f"Field '{field}' is not defined in configuration.")

    field_type = field_config.type
    if field_type == "calculated":
        # Get the fully qualified data column table.column specifier
        data_column_spec = f"{data_column.table.name}.{data_column.key}"

        # Substitute keys in the expression with their JSONB paths
        expression: str = field_config.expression
        for key, prop in properties_config.items():
            if prop.type != 'calculated':
                json_key = prop.json_key
                expression = expression.replace(
                    f"{{{key}}}",
                    f"(CAST({data_column_spec}::jsonb->>'{json_key}' AS FLOAT))"
                )
        # Substitute parameters if provided
        if params:
            for param_key, param_value in params.items():
                expression = expression.replace(f"{{{param_key}}}", str(param_value))
        return literal_column(expression)
    elif field_type == "float":
        json_key = field_config.json_key
        return cast(data_column.op('->>')(json_key), Float)
    else:
        raise ValueError(f"Unsupported field type '{field_type}'.")