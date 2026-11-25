from typing import Annotated

from fastapi import Depends
from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import Session


# pass empty connection string to use PG* environment variables (see https://www.postgresql.org/docs/current/libpq-envars.html)
engine = create_engine("postgresql+psycopg2://", future=True, pool_pre_ping=True)


def get_session():
    with Session(engine) as session:
        yield session


SessionDep = Annotated[Session, Depends(get_session)]

Base = declarative_base()
