from sqlalchemy import ForeignKey, Integer, Column, String, JSON, Float
from sqlalchemy.orm import relationship
from geoalchemy2 import Geometry

from .database import Base


class Feature(Base):
    __tablename__ = "features"

    id = Column(String, primary_key=True)
    layer = Column(String, index=True, nullable=False)
    sublayer = Column(String)
    properties = Column(JSON, nullable=False)
    geom = Column(Geometry, nullable=False)

    damages = relationship("Damage")


class Damage(Base):
    __tablename__ = "damages"

    feature_id = Column("id", String, ForeignKey(Feature.id), primary_key=True)
    hazard = Column(String, nullable=False, primary_key=True)
    rcp = Column(String, nullable=False, primary_key=True)
    epoch = Column(String, nullable=False, primary_key=True)

    EAD_undefended_min = Column(Float)
    EAD_undefended_mean = Column(Float)
    EAD_undefended_max = Column(Float)

    EAEL_undefended_min = Column(Float)
    EAEL_undefended_mean = Column(Float)
    EAEL_undefended_max = Column(Float)
