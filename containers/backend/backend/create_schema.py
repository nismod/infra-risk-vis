"""
SQLAlchemy Schema Setup Main
"""

from db import models
from db.database import engine

if __name__ == "__main__":
    models.Base.metadata.drop_all(bind=engine)
    models.Base.metadata.create_all(bind=engine)
