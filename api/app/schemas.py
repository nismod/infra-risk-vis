from pydantic import BaseModel


class FeatureBase(BaseModel):
    id: str
    layer: str
    sublayer: str
    properties: dict


class Feature(FeatureBase):
    class Config:
        orm_mode = True
