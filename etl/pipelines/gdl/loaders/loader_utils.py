import json
import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker


# Not needed in prod when gri_metrics package (defined in containers/metrics/pyproject.toml) is installed in the working python environment
# Keeping here because useful when setting up a dev env
def insert_paths():
    import sys
    import os
    import inspect

    current_dir = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe()))
    )
    parent_dir = os.path.dirname(os.path.dirname(current_dir))
    grand_dir = os.path.dirname(os.path.dirname(parent_dir))
    etl_dir = os.path.dirname(os.path.dirname(grand_dir))
    container_dir = os.path.dirname(os.path.dirname(etl_dir))
    sys.path.insert(0, parent_dir)
    sys.path.insert(0, grand_dir)
    sys.path.insert(0, etl_dir)
    sys.path.insert(0, container_dir)


# insert_paths()


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
