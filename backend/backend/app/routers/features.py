from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from backend.app import schemas
from backend.app.dependencies import get_db
from backend.db import models


router = APIRouter(tags=["features"])


@router.get("/{feature_id}", response_model=schemas.Feature)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    return feature
