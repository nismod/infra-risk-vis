from sqlalchemy.orm import declarative_base
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
import sqlalchemy as sa


import os
from dotenv import load_dotenv

load_dotenv()


# Shared with ETL
def get_db_uri() -> sa.engine.URL:
    return os.getenv("METRICS_DB_URL")


engine = create_engine(get_db_uri(), pool_pre_ping=True)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
