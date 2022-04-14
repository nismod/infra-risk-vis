from sqlalchemy import ForeignKey, Integer, Column, String, JSON, Float
from sqlalchemy.orm import relationship
from geoalchemy2 import Geometry

from .database import Base


class Feature(Base):
    __tablename__ = "features"
    id = Column(Integer, primary_key=True)
    string_id = Column(String, nullable=False)
    layer = Column(String, index=True, nullable=False)
    sublayer = Column(String)
    properties = Column(JSON, nullable=False)
    geom = Column(Geometry("GEOMETRY", srid=4326), nullable=False)

    damages = relationship("Damage")


class Damage(Base):
    __tablename__ = "damages"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True)

    hazard = Column(String, nullable=False, primary_key=True)
    rcp = Column(String, nullable=False, primary_key=True)
    epoch = Column(String, nullable=False, primary_key=True)

    # Damage type (direct, indirect, combined)
    damage_type = Column(String, nullable=False, primary_key=True)

    protection_standard = Column(Integer, nullable=False, primary_key=True)

    min = Column(Float)
    mean = Column(Float)
    max = Column(Float)
