import json
import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker
from containers.metrics.api.database.database import get_db_uri


def load_json(fpath: str):
    with open(fpath, "r") as file:
        data = json.load(file)

    return data


def init_db_session():
    db_uri = get_db_uri()
    engine = sa.create_engine(db_uri, pool_pre_ping=True)
    Session = sessionmaker(autocommit=False, autoflush=False, bind=engine)

    return engine, Session
