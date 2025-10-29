"""DB Models to Support GDL API"""

from sqlalchemy import Column, String, SmallInteger, Float, ForeignKey
from sqlalchemy.orm import relationship, declared_attr
from geoalchemy2 import Geometry
from .database import Base


class Region(Base):
    """Generic regional boundaries table"""

    __tablename__ = "region"

    source_code = Column(String, index=True)
    region_name = Column(String, nullable=False)
    level = Column(String, nullable=False)
    geometry = Column(Geometry("GEOMETRY", srid=4326), nullable=False)
    iso_code = Column(String, ForeignKey("iso_country.iso_code"), nullable=False)
    iso_country = relationship("IsoCountry", back_populates="region")


class IsoCountry(Base):
    """ISO country codes table"""

    __tablename__ = "iso_country_codes"

    iso_code = Column(String, primary_key=True, index=True)
    country_name = Column(String, nullable=False)
    continent = Column(String, nullable=False)
    regions = relationship("Region", back_populates="iso_country")


class Metric(Base):
    """Region metrics"""

    __tablename__ = "region_metric"

    source_code = Column(String, ForeignKey("region.source_code"), index=True)
    region = relationship("MetricRegion")

    year = Column(SmallInteger, index=True)

    development = Column(Float, nullable=True)
    education = Column(Float, nullable=True)
    healthcare = Column(Float, nullable=True)
    income = Column(Float, nullable=True)
