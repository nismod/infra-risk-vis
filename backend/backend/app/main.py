from fastapi import FastAPI, Depends, Path, Query
from pydantic import Json
from sqlalchemy.orm import Session

from backend.db import models
from backend.db.database import engine

from . import schemas
from .dependencies import get_db
from .routers import attributes


app = FastAPI()


@app.get("/features/{feature_id}", response_model=schemas.Feature)
def read_feature(feature_id: int, db: Session = Depends(get_db)):
    feature = db.query(models.Feature).filter(models.Feature.id == feature_id).one()
    return feature


app.include_router(attributes.router, prefix="/attributes")
