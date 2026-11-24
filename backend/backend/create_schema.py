"""
SQLAlchemy Schema Setup Main
"""

from backend.db import models
from backend.db.database import engine

if __name__ == "__main__":
    models.Base.metadata.drop_all(bind=engine)
    models.Base.metadata.create_all(bind=engine)
