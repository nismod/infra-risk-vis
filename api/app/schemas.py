from pydantic import BaseModel


class FeatureBase(BaseModel):
    id: str
    layer: str
    sublayer: str | None
    properties: dict


class DamageBase(BaseModel):
    hazard: str
    rcp: str
    epoch: str
    EAD_undefended_min: float
    EAD_undefended_mean: float
    EAD_undefended_max: float
    EAEL_undefended_min: float
    EAEL_undefended_mean: float
    EAEL_undefended_max: float


class Damage(DamageBase):
    class Config:
        orm_mode = True


class Feature(FeatureBase):
    damages: list[Damage] = []

    class Config:
        orm_mode = True
