from enum import Enum
from typing import Any, Generic, Literal, TypeVar
from pydantic import BaseModel, conint, root_validator
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


class DataDimensions(BaseModel):
    pass


class DataVariables(BaseModel):
    pass


class DataParameters(BaseModel):
    pass


# Expected Damages
class ExpectedDamagesDimensions(DataDimensions):
    hazard: str
    rcp: str
    epoch: str
    protection_standard: int


class ExpectedDamagesVariables(DataVariables):
    ead_amin: float
    ead_mean: float
    ead_amax: float
    eael_amin: float
    eael_mean: float
    eael_amax: float


class ExpectedDamage(ExpectedDamagesDimensions, ExpectedDamagesVariables):
    class Config:
        orm_mode = True


# Return Period Damages
class ReturnPeriodDamagesDimensions(DataDimensions):
    hazard: str
    rcp: str
    epoch: str
    rp: int


class ReturnPeriodDamagesVariables(DataVariables):
    exposure: float
    damage_amin: float
    damage_mean: float
    damage_amax: float
    loss_amin: float
    loss_mean: float
    loss_amax: float


class ReturnPeriodDamage(ReturnPeriodDamagesDimensions, ReturnPeriodDamagesVariables):
    class Config:
        orm_mode = True


# NPV Damages
class NPVDamagesDimensions(DataDimensions):
    hazard: str
    rcp: str


class NPVDamagesVariables(DataVariables):
    ead_amin: float
    ead_mean: float
    ead_amax: float
    eael_amin: float
    eael_mean: float
    eael_amax: float


class NPVDamage(NPVDamagesDimensions, NPVDamagesVariables):
    class Config:
        orm_mode = True


# Adaptation Options


class AdaptationDimensions(DataDimensions):
    hazard: str
    rcp: str
    adaptation_name: str
    adaptation_protection_level: float


class AdaptationVariables(DataVariables):
    adaptation_cost: float

    avoided_ead_amin: float
    avoided_ead_mean: float
    avoided_ead_amax: float
    avoided_eael_amin: float
    avoided_eael_mean: float
    avoided_eael_amax: float


class AdaptationCostBenefitRatioParameters(DataParameters):
    eael_days: conint(ge=1, le=30)


class Adaptation(AdaptationDimensions, AdaptationVariables):
    class Config:
        orm_mode = True


# Features
class FeatureOutBase(FeatureBase):
    class Config:
        orm_mode = True


class FeatureOut(FeatureOutBase):
    damages_expected: list[ExpectedDamage] = []
    damages_return_period: list[ReturnPeriodDamage] = []
    damages_npv: list[NPVDamage] = []
    adaptations: list[Adaptation] = []


# Features Sorted Lists


class LayerSpec(BaseModel):
    layer_name: str | None
    sector: str | None
    subsector: str | None
    asset_type: str | None


SortFieldT = TypeVar("SortFieldT")


class FeatureListItemOut(GenericModel, Generic[SortFieldT]):
    id: int
    string_id: str
    layer: str
    bbox_wkt: str
    value: SortFieldT

    class Config:
        orm_mode = True


# Feature Attributes Lookups

AttributeT = TypeVar("AttributeT")

AttributeLookup = dict[int, AttributeT]
