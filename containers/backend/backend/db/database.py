from sqlalchemy import create_engine
from sqlalchemy.orm import declarative_base
from sqlalchemy.orm import sessionmaker


# pass empty connection string to use PG* environment variables (see https://www.postgresql.org/docs/current/libpq-envars.html)
engine = create_engine("postgresql+psycopg2://", future=True, pool_pre_ping=True)

SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
