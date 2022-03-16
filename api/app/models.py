from sqlalchemy import Integer, Column, String, JSON
from geoalchemy2 import Geometry

from .database import Base


class Feature(Base):
    __tablename__ = "features"

    id = Column(String, primary_key=True, index=True)
    layer = Column(String, index=True)
    sublayer = Column(String)
    properties = Column(JSON)
    geom = Column(Geometry)
