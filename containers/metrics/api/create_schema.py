"""
SQLAlchemy Schema Setup Main
"""

import database.gdl as gdl
from api.database.database import engine

if __name__ == "__main__":
    gdl.Base.metadata.drop_all(bind=engine)
    gdl.Base.metadata.create_all(bind=engine)
