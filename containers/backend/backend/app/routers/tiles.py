"""
Tile Service wrapping Terracotta Python API
"""
from sys import getsizeof
from typing import Any, BinaryIO, List, Tuple, Union
import ast
import inspect

from fastapi import APIRouter, Depends, HTTPException, Response, Header
from fastapi.logger import logger
from fastapi.responses import FileResponse
import sqlalchemy
from starlette.responses import StreamingResponse
from mercantile import tile
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import Session
from geoalchemy2 import functions


from app import schemas
from app.dependencies import get_db
from db import models
from app.internal.helpers import build_driver_path, handle_exception
from app.exceptions import SourceDBDoesNotExistException, DomainAlreadyExistsException
from config import API_TOKEN


router = APIRouter(tags=["tiles"])


def _parse_keys(keys: str) -> List:
    """
    Parse tiles URL key str
    """
    return [key for key in keys.split("/") if key]


def _domain_from_keys(keys: str) -> str:
    """
    Retrieve domain from keys path
    """
    return _parse_keys(keys)[0]


def _get_singleband_image(
    database: str, keys: str, tile_xyz: Tuple[int, int, int] = None, options: dict = {}
) -> BinaryIO:

    """
    Generate a Singleband Tile

    ::args database str DB under-which the requested data has been loaded
    """
    from app.internal.tiles.singleband import singleband

    parsed_keys = _parse_keys(keys)

    if options.get("colormap", "") == "explicit":
        options["colormap"] = options.pop("explicit_color_map")

    # Collect TC Driver path for MySQL
    driver_path = build_driver_path(database)

    logger.debug(
        "driver path: %s, parsed_keys: %s, tile_xyz: %s, options: %s",
        driver_path,
        parsed_keys,
        tile_xyz,
        options,
    )

    return singleband(driver_path, parsed_keys, tile_xyz=tile_xyz, **options)


def _source_db_exists(db: Session, source_db: str) -> bool:
    """
    Check whether the given source_db exists in the meta store
    """
    res = (
        db.query(models.RasterTileSource)
        .filter(models.RasterTileSource.source_db == source_db)
        .all()
    )
    if res:
        return True
    return False


def _domain_exists(db: Session, domain: str) -> bool:
    """
    Check whether the given domain exists in the meta store
    """
    res = (
        db.query(models.RasterTileSource)
        .filter(models.RasterTileSource.domain == domain)
        .all()
    )
    if res:
        return True
    return False


def _tile_db_from_keys(db: Session, keys: List) -> str:
    """
    Query the name of the mysql database within-which the tiles reside
        using the first value of the path keys (which links to hazard-type in the UI)
        See: frontend/src/config/hazards/domains.ts
    """
    domain = _domain_from_keys(keys)
    # Should only return one entry
    res = (
        db.query(models.RasterTileSource)
        .filter(models.RasterTileSource.domain == domain)
        .one()
    )
    return res.source_db


async def verify_token(x_token: str = Header()):
    print(x_token, API_TOKEN)
    if x_token != API_TOKEN:
        raise HTTPException(status_code=401, detail="")


@router.get("/sources", response_model=List[schemas.TileSourceMeta])
def get_all_tile_source_meta(
    db: Session = Depends(get_db),
) -> List[schemas.TileSourceMeta]:
    """
    Retrieve metadata about all the tile sources available
    """
    logger.debug(
        "performing %s",
        inspect.stack()[0][3],
    )
    try:
        res = db.query(models.RasterTileSource).all()
        return [schemas.TileSourceMeta.from_orm(row) for row in res]
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.get("/sources/{source_id}", response_model=schemas.TileSourceMeta)
def get_tile_source_meta(
    source_id: int,
    db: Session = Depends(get_db),
) -> List[schemas.TileSourceMeta]:
    """
    Retrieve metadata about a single tile source
    """
    logger.debug(
        "performing %s",
        inspect.stack()[0][3],
    )
    try:
        res = (
            db.query(models.RasterTileSource)
            .filter(models.RasterTileSource.id == source_id)
            .one()
        )
        return res
    except NoResultFound:
        raise HTTPException(status_code=404)
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.post("/sources", dependencies=[Depends(verify_token)])
def insert_source_meta(
    source_meta: schemas.TileSourceMeta, db: Session = Depends(get_db)
) -> Any:
    """
    Ingest Tile Source Meta
    """
    logger.debug("performing %s, with %s", inspect.stack()[0][3], source_meta)
    try:
        # Check if the domain already exists
        if _domain_exists(db, source_meta.domain):
            raise DomainAlreadyExistsException()
        values = source_meta.dict()
        values.pop("id")
        stmt = sqlalchemy.insert(models.RasterTileSource).values(values)
        res = db.execute(stmt)
        db.commit()
        logger.debug("Insert Meta result: %s", res)
        return Response(status_code=200)
    except DomainAlreadyExistsException as err:
        handle_exception(logger, err)
        raise HTTPException(
            status_code=400, detail="domain already exists in another database"
        )
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.delete("/sources/{source_id:int}", dependencies=[Depends(verify_token)])
def delete_source_meta(source_id: int, db: Session = Depends(get_db)) -> Any:
    """
    Delete Tile Source Meta
    """
    logger.debug("performing %s, with %s", inspect.stack()[0][3], source_id)
    try:
        stmt = sqlalchemy.delete(models.RasterTileSource).filter(
            models.RasterTileSource.id == source_id
        )
        res = db.execute(stmt)
        db.commit()
        logger.debug("Delete Meta result: %s", res)
        return Response(status_code=200)
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.get(
    "/{keys:path}/{tile_z:int}/{tile_x:int}/{tile_y:int}.png",
    response_class=StreamingResponse,
)
async def get_tile(
    keys: str,
    tile_z: int,
    tile_x: int,
    tile_y: int,
    colormap: Union[str, None] = None,
    stretch_range: Union[str, None] = None,
    db: Session = Depends(get_db),
):
    """
    This route does not work in a FastAPI thread pool environment (i.e. when not async)
    """
    logger.debug(
        "tile path %s, colormap: %s, stretch_range: %s",
        keys,
        colormap,
        ast.literal_eval(stretch_range),
    )
    try:
        # Collect the source tiles DB based on the first value of the path (the domain or type)
        source_db = _tile_db_from_keys(db, keys)
        logger.debug("source DB for tile path: %s", source_db)
    except NoResultFound:
        domain = _domain_from_keys(keys)
        raise HTTPException(
            status_code=404,
            detail=f"no source database for the given domain {domain} could be found - was the source metadata loaded to PG?",
        )
    try:
        # Check the keys are appropriate for the given type
        options = {}
        if colormap:
            options["colormap"] = colormap
        if stretch_range:
            options["stretch_range"] = ast.literal_eval(stretch_range)

        # Check the database exists
        if not _source_db_exists(db, source_db):
            raise SourceDBDoesNotExistException()

        # Generate the tile
        image = _get_singleband_image(
            source_db, keys, [tile_x, tile_y, tile_z], options
        )
        logger.debug(
            "tile image of size returned: %s, %s", getsizeof(image), type(image)
        )

        # Return the tile as stream
        return StreamingResponse(image, media_type="image/png")
    except SourceDBDoesNotExistException as err:
        handle_exception(logger, err)
        raise HTTPException(
            status_code=400,
            detail=f"source database {source_db} does not exist in tiles metastore",
        )
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)