from enum import Enum
from typing import Generic, Optional, TypeVar, Union
from typing_extensions import Annotated
from pydantic import BaseModel, ConfigDict, Field, field_validator


class FeatureBase(BaseModel):
    id: int
    string_id: str
    layer: str
    sublayer: Optional[str] = None
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
    epoch: Union[str, int]
    protection_standard: int


class ExpectedDamagesVariables(DataVariables):
    ead_amin: float
    ead_mean: float
    ead_amax: float
    eael_amin: float
    eael_mean: float
    eael_amax: float


class ExpectedDamage(ExpectedDamagesDimensions, ExpectedDamagesVariables):
    model_config = ConfigDict(from_attributes=True)


# Return Period Damages
class ReturnPeriodDamagesDimensions(DataDimensions):
    hazard: str
    rcp: str
    epoch: Union[str, int]
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
    model_config = ConfigDict(from_attributes=True)


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
    model_config = ConfigDict(from_attributes=True)


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
    eael_days: Annotated[int, Field(ge=1, le=30)]

    @field_validator("eael_days")
    def fix_eael_days(cls, eael_days: int) -> float:
        """
        The data for `AdaptationCostBenefit.avoided_eael_mean` is erroneous and
        should be modified in the meantime. This validator adds a fudge factor
        to account for the multiplicative error in the data.

        TODO: fix the data so we don't need this function.
        """
        return eael_days / 15


class Adaptation(AdaptationDimensions):
    properties: dict
    model_config = ConfigDict(from_attributes=True)


# Features
class FeatureOutBase(FeatureBase):
    model_config = ConfigDict(from_attributes=True)


class FeatureOut(FeatureOutBase):
    damages_expected: list[ExpectedDamage] = []
    damages_return_period: list[ReturnPeriodDamage] = []
    damages_npv: list[NPVDamage] = []
    adaptation: list[Adaptation] = []


# Features Sorted Lists


class LayerSpec(BaseModel):
    layer_name: str | None
    sector: str | None
    subsector: str | None
    asset_type: str | None


SortFieldT = TypeVar("SortFieldT")


class FeatureListItemOut(BaseModel, Generic[SortFieldT]):
    id: int
    string_id: str
    layer: str
    bbox_wkt: str
    value: SortFieldT

    model_config = ConfigDict(from_attributes=True)


# Feature Attributes Lookups

AttributeT = TypeVar("AttributeT")

AttributeLookup = dict[int, AttributeT]


# Tile Server metadata
class TileSourceMeta(BaseModel):
    id: Optional[int] = None
    domain: str
    name: str
    group: str
    description: str
    license: str
    keys: list[str]

    model_config = ConfigDict(from_attributes=True)


class TileSourceDomains(BaseModel):
    domains: list[dict[str, str]]


class ColorMapOptions(BaseModel):
    stretch_range: list[int]
    colormap: str
    num_values: Optional[int] = 255


class ColorMapEntry(BaseModel):
    value: float
    rgba: list[int]


class ColorMap(BaseModel):
    colormap: list[ColorMapEntry]


# Pixel Driller Schemas
class PixelDrillerLayer(BaseModel):
    domain: str
    keys: dict[str, Union[str, int, float]]


class PixelDrillerResult(BaseModel):
    value: Optional[float] = None
    layer: PixelDrillerLayer


class PixelDrillerResponse(BaseModel):
    point: dict[str, float]
    results: list[PixelDrillerResult]
