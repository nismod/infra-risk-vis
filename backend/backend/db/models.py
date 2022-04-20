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

    damages_return_period = relationship("ReturnPeriodDamage")
    damages_expected = relationship("ExpectedDamage")
    damages_npv = relationship("NPVDamage")
    adaptation = relationship("AdaptationCostBenefit")


class ReturnPeriodDamage(Base):
    __tablename__ = "damages_rp"

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True)

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

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True)

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

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True)

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

    feature_id = Column(Integer, ForeignKey(Feature.id), primary_key=True)

    hazard = Column(String(8), nullable=False, primary_key=True)
    rcp = Column(String(8), nullable=False, primary_key=True)

    adaptation_name = Column(String, nullable=False, primary_key=True)
    adaptation_protection_level = Column(Float, nullable=False, primary_key=True)

    adaptation_cost = Column(Float)

    avoided_ead_amin = Column(Float)
    avoided_ead_mean = Column(Float)
    avoided_ead_amax = Column(Float)

    avoided_eael_amin = Column(Float)
    avoided_eael_mean = Column(Float)
    avoided_eael_amax = Column(Float)
