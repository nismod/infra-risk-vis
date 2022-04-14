from enum import Enum
from typing import Any, Generic, Literal, TypeVar
from pydantic import BaseModel
from pydantic.generics import GenericModel


class FeatureBase(BaseModel):
    id: int
    string_id: str
    layer: str
    sublayer: str | None
    properties: dict


class DamageType(str, Enum):
    direct = "direct"
    indirect = "indirect"
    combined = "combined"


class DamageParams(BaseModel):
    hazard: str
    rcp: str
    epoch: str
    damage_type: DamageType
    protection_standard: int


class DamageBase(DamageParams):
    min: float
    mean: float
    max: float


class Damage(DamageBase):
    class Config:
        orm_mode = True


class Feature(FeatureBase):
    damages: list[Damage] = []

    class Config:
        orm_mode = True


AttributeT = TypeVar("AttributeT")

AttributeLookup = dict[int, AttributeT]
