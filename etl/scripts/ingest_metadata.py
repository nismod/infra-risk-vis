import json

from sqlalchemy import insert, delete
from sqlalchemy.orm import Session
from backend.db.models import RasterTileSource
from backend.db.database import SessionLocal


def create_source_meta(meta: dict, db: Session):
    stmt = delete(RasterTileSource).where(RasterTileSource.domain == meta["domain"])
    db.execute(stmt)
    stmt = insert(RasterTileSource).values(meta)
    db.execute(stmt)
    db.commit()


if __name__ == "__main__":
    try:
        path = snakemake.input.metadata
        flag = snakemake.output.flag
    except NameError:
        assert False, "Must be run from snakemake"

    with open(path) as fh:
        metadata = json.load(fh)
    db = SessionLocal()
    create_source_meta(metadata, db)

    with open(flag, "w") as fh:
        fh.write("Done")
