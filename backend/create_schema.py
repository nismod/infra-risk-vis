from backend.db import models
from backend.db.database import engine

models.Base.metadata.create_all(bind=engine)
