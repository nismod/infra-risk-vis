"""
SQL Alchemy Models for Backend
"""

from sqlalchemy import ForeignKey, Integer, Column, String, JSON, Float
from sqlalchemy.dialects.postgresql import JSONB
from sqlalchemy.orm import relationship
from geoalchemy2 import Geometry

from .database import Base


class FeatureLayer(Base):
    __tablename__ = "feature_layers"
    layer_name = Column(String, primary_key=True)
    sector = Column(String, nullable=False)
    subsector = Column(String, nullable=False)
    asset_type = Column(String, nullable=False)


class Feature(Base):
    __tablename__ = "features"
    id = Column(Integer, primary_key=True)
    string_id = Column(String, nullable=False)
    layer = Column(
        String, ForeignKey(FeatureLayer.layer_name), index=True, nullable=False
    )
    properties = Column(JSONB, nullable=False)
    geom = Column(Geometry("GEOMETRY", srid=4326), nullable=False)

    layer_info = relationship("FeatureLayer")

    damages_return_period = relationship("ReturnPeriodDamage")
    damages_expected = relationship("ExpectedDamage")
    damages_npv = relationship("NPVDamage")
    adaptation = relationship("AdaptationCostBenefit")


class ReturnPeriodDamage(Base):
    __tablename__ = "damages_rp"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True, index=True)

    hazard = Column(String(8), nullable=False, primary_key=True)
    rcp = Column(String(8), nullable=False, primary_key=True)
    epoch = Column(Integer, nullable=False, primary_key=True)
    rp = Column(Integer, nullable=False, primary_key=True)

    exposure = Column(Float)

    damage_amin = Column(Float)
    damage_mean = Column(Float)
    damage_amax = Column(Float)

    loss_amin = Column(Float)
    loss_mean = Column(Float)
    loss_amax = Column(Float)


class ExpectedDamage(Base):
    __tablename__ = "damages_expected"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True, index=True)

    hazard = Column(String(8), nullable=False, primary_key=True)
    rcp = Column(String(8), nullable=False, primary_key=True)
    epoch = Column(Integer, nullable=False, primary_key=True)

    protection_standard = Column(Integer, nullable=False, primary_key=True)

    ead_amin = Column(Float)
    ead_mean = Column(Float)
    ead_amax = Column(Float)

    eael_amin = Column(Float)
    eael_mean = Column(Float)
    eael_amax = Column(Float)


class NPVDamage(Base):
    __tablename__ = "damages_npv"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True, index=True)

    hazard = Column(String(8), nullable=False, primary_key=True)
    rcp = Column(String(8), nullable=False, primary_key=True)

    ead_amin = Column(Float)
    ead_mean = Column(Float)
    ead_amax = Column(Float)

    eael_amin = Column(Float)
    eael_mean = Column(Float)
    eael_amax = Column(Float)


class AdaptationCostBenefit(Base):
    __tablename__ = "adaptation_cost_benefit"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True, index=True)

    hazard = Column(String(8), nullable=False, primary_key=True)
    rcp = Column(String(8), nullable=False, primary_key=True)

    adaptation_name = Column(String, nullable=False, primary_key=True)
    adaptation_protection_level = Column(Float, nullable=False, primary_key=True)

    properties = Column(JSONB, nullable=False)



class RasterTileSource(Base):
    """
    Can be multiple RasterTileSources to a single source mysql databae
    """

    __tablename__ = "raster_tile_sources"
    id = Column(Integer, primary_key=True, autoincrement=True)
    # Domain of the TileSource (this is used to check which database a front-end call goes to)
    # Domains can only reside in a single database - we dont allow duplicates even between databases
    domain = Column(String, nullable=False, unique=True)
    name = Column(String, nullable=False)
    group = Column(String, nullable=False)  # Hazard, Risk, Exposure, Adaptation
    description = Column(String, nullable=True)
    license = Column(String, nullable=True)
    keys = Column(JSONB)  # JSON list of terracotta/URL keys
