from fastapi import Query
from backend.app import schemas
from backend.db.database import SessionLocal


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def get_damage_params(
    hazard: str,
    rcp: str,
    epoch: str,
    damage_type: schemas.DamageType,
    protection_standard: int,
):
    args = locals()  # https://stackoverflow.com/a/2521937/1478817
    return schemas.DamageParams(**args)
