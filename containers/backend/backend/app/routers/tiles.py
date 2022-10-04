"""
Tile Service wrapping Terracotta Python API
"""
from sys import getsizeof
from typing import Any, BinaryIO, List, Tuple, Union
import ast
import inspect

from fastapi import APIRouter, Depends, HTTPException, Response
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
from app.exceptions import SourceDBAlreadyExistsException

router = APIRouter(tags=["tiles"])


def _get_singleband_image(
    database: str, keys: str, tile_xyz: Tuple[int, int, int] = None, options: dict = {}
) -> BinaryIO:

    """
    Generate a Singleband Tile

    ::args database str DB under-which the requested data has been loaded
    """
    from app.internal.tiles.singleband import singleband

    parsed_keys = [key for key in keys.split("/") if key]

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


@router.post("/sources")
def insert_source_meta(
    source_meta: schemas.TileSourceMeta, db: Session = Depends(get_db)
) -> Any:
    """
    Ingest Tile Source Meta
    """
    logger.debug("performing %s, with %s", inspect.stack()[0][3], source_meta)
    try:
        res = (
            db.query(models.RasterTileSource)
            .filter(models.RasterTileSource.source_db == source_meta.source_db)
            .all()
        )
        if res:
            raise SourceDBAlreadyExistsException()
        values = source_meta.dict()
        values.pop("id")
        stmt = sqlalchemy.insert(models.RasterTileSource).values(values)
        res = db.execute(stmt)
        db.commit()
        logger.debug("Insert Meta result: %s", res)
        return Response(status_code=200)
    except SourceDBAlreadyExistsException as err:
        handle_exception(logger, err)
        raise HTTPException(
            status_code=400,
            detail=f"the source tiles database {source_meta.source_db} already exists",
        )
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.delete("/sources/{source_id:int}")
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
    "/{source_db:str}/{keys:path}/{tile_z:int}/{tile_x:int}/{tile_y:int}.png",
    response_class=StreamingResponse,
)
async def get_tile(
    source_db: str,
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
        "tile db_source: %s, tile path %s, colormap: %s, stretch_range: %s",
        source_db,
        keys,
        colormap,
        ast.literal_eval(stretch_range),
    )

    # Check the keys are appropriate for the given type
    options = {}
    if colormap:
        options["colormap"] = colormap
    if stretch_range:
        options["stretch_range"] = ast.literal_eval(stretch_range)

    # Generate the tile
    image = _get_singleband_image(source_db, keys, [tile_x, tile_y, tile_z], options)
    logger.debug("tile image of size returned: %s, %s", getsizeof(image), type(image))

    # Return the tile as stream
    return StreamingResponse(image, media_type="image/png")
