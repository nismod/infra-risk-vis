"""
Tile Service wrapping Terracotta Python API
"""
from collections import OrderedDict
from sys import getsizeof
from typing import Any, BinaryIO, List, Tuple, Union
import ast
import inspect
import itertools

from fastapi import APIRouter, Depends, HTTPException, Response, Header
from fastapi.logger import logger
import sqlalchemy
from starlette.responses import StreamingResponse
from sqlalchemy.exc import NoResultFound
from sqlalchemy.orm import Session
from geoalchemy2 import functions


from app import schemas
from app.dependencies import get_db
from db import models
from app.internal.tiles.singleband import all_datasets, database_keys
from app.internal.helpers import build_driver_path, handle_exception
from app.exceptions import (
    SourceDBDoesNotExistException,
    DomainAlreadyExistsException,
    MissingExplicitColourMapException,
)
from config import API_TOKEN, DOMAIN_TO_DB_MAP, CATEGORICAL_COLOR_MAPS


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


def _get_tiledb_keys(database: str) -> List[str]:
    """Retrieve keys and ordering for given DB"""
    driver_path = build_driver_path(database)
    db_keys = database_keys(driver_path)
    return list(db_keys.keys())


def _get_singleband_image(
    database: str, keys: str, tile_xyz: Tuple[int, int, int] = None, options: dict = {}
) -> BinaryIO:

    """
    Generate a Singleband Tile

    ::args database str DB under-which the requested data has been loaded
    """
    from app.internal.tiles.singleband import singleband

    parsed_keys = _parse_keys(keys)

    # Collect TC Driver path for MySQL
    driver_path = build_driver_path(database)

    logger.debug(
        "parsed_keys: %s, tile_xyz: %s, options: %s",
        parsed_keys,
        tile_xyz,
        options,
    )

    return singleband(driver_path, parsed_keys, tile_xyz=tile_xyz, **options)


def _source_db_exists(source_db: str) -> bool:
    """
    Check whether the given source_db exists in the meta store
    """
    return source_db in DOMAIN_TO_DB_MAP.values()


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


def _tile_db_from_keys(keys: List) -> str:
    """
    Query the name of the mysql database within-which the tiles reside
        using the first value of the path keys (which links to hazard-type in the UI)
        See: frontend/src/config/hazards/domains.ts
    """
    domain = _domain_from_keys(keys)
    # Should only return one entry
    return DOMAIN_TO_DB_MAP[domain]


def _source_options(source_db: str, domain: str = None) -> List[dict]:
    """
    Gather all URL key combinations available in the given source

    ::param source_db str The name of the source MySQL Database in-which the tiles reside
    ::kwarg domain str Source options will be optionally filtered to only include this 'type' (the first index in the key-tuple)
        domain in the source meta should always map to the type - which is the first key in the key tuple
        e.g. 'fluvial in th example tile k:v pairs below'

    ('fluvial', '50', '8x5', '2080', 'NorESM1-M'):
        '/data/aqueduct/inunriver_rcp8p5_00000NorESM1-M_2080_rp00050.tif',
    ('fluvial', '50', 'baseline', 'present', 'WATCH'):
        '/data/aqueduct/inunriver_historical_000000000WATCH_1980_rp00050.tif',
    ('fluvial', '500', '4x5', '2030', 'GFDL-ESM2M'):
        '/data/aqueduct/inunriver_rcp4p5_0000GFDL-ESM2M_2030_rp00500.tif',
    ('fluvial', '500', '4x5', '2030', 'HadGEM2-ES'):
        '/data/aqueduct/inunriver_rcp4p5_0000HadGEM2-ES_2030_rp00500.tif'
    """

    driver_path = build_driver_path(source_db)

    # query terracotta for all datasets in the source
    # dict of tuples of variable values
    # e.g. ("2020", "10", "constant") pointing to raster path strs
    datasets: dict[tuple, str] = all_datasets(driver_path)

    # get the key names for these data
    # e.g. {'epoch': '', 'rp': '', 'ssp': ''}
    keys: OrderedDict[str, str] = database_keys(driver_path)

    # now generate the output mapping
    # take the variable values from datasets and create dicts naming them as such
    # e.g. [{"epoch": "2020", "rp": "10", "ssp": "constant"}, ...]
    source_options: list[dict[str, str]] = [
        dict(zip(keys, _values)) for _values in datasets.keys()
    ]

    logger.debug(f"{source_db=} {domain=} {driver_path=} {datasets=} {keys=} {source_options=}")

    # optionally filter to a domain (type)
    if domain is not None:
        source_options = [item for item in source_options if item["type"] == domain]

    return source_options


async def verify_token(x_token: str = Header()):
    if x_token != API_TOKEN:
        raise HTTPException(status_code=401, detail="")


@router.get("/sources", response_model=List[schemas.TileSourceMeta])
async def get_all_tile_source_meta(
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
        all_meta = []
        for row in res:
            url_keys = _get_tiledb_keys(row.source_db)
            meta = schemas.TileSourceMeta.from_orm(row)
            meta.url_keys = url_keys
            all_meta.append(meta)
        return all_meta
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.get("/sources/{source_id}", response_model=schemas.TileSourceMeta)
async def get_tile_source_meta(
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
        url_keys = _get_tiledb_keys(res.source_db)
        meta = schemas.TileSourceMeta.from_orm(res)
        meta.url_keys = url_keys
        return meta
    except NoResultFound:
        raise HTTPException(status_code=404)
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.get("/sources/{source_id}/domains", response_model=schemas.TileSourceDomains)
async def get_tile_source_domains(
    source_id: int,
    db: Session = Depends(get_db),
) -> schemas.TileSourceDomains:
    """
    Retrieve all combinations available for the source domain
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
        domains = _source_options(
            res.source_db, domain=res.domain if res.domain else None
        )
        meta = schemas.TileSourceDomains(domains=domains)
        logger.debug(f"{source_id=} {res.source_db=} {res.domain=} {domains=} {meta=}")
        return meta
    except NoResultFound:
        raise HTTPException(status_code=404)
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)


@router.post("/sources", dependencies=[Depends(verify_token)])
async def insert_source_meta(
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
        values.pop("url_keys")
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
async def delete_source_meta(source_id: int, db: Session = Depends(get_db)) -> Any:
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
    explicit_color_map: Union[str, None] = None,
    db: Session = Depends(get_db),
):
    """
    Serves XYZ Raster Tiles with the given colormap / stretch range or explicit color map for categorical data.

    ::param keys str A string containing the url-encoded keys which address the required raster.

        This string is constructed as-per Terracotta URL requests, e.g: `aqueduct/gsm1/1980/baseline`

        Information about datasets available, their associated keys and order can be found using the `/tiles/sources/21/domains` endpoint.

    ::param tile_z int Tile Z address
    ::param tile_x int Tile X address
    ::param tile_y int Tile Y address

    ::kwarg colormap str A string representing the colormap to be used to render the tile.  Colormaps can be access separately through the `/colormap` endpoint

        e.g. `colormap=reds`

    ::kwarg stretch_range iterable The range over-which to stretch the pixel values

        e.g. `stretch_range=[0,10]`

    ::kwarg explicit_color_map str A categorical colormap `{pixel_value: (R,G,B,A)}` to be used with a given categorical data-source.

        __NOTE__: `colormap` arg must be set to "explicit" in order to use `explicit_color_map`

        .e.g colormap=explicit&explicit_color_map="{\"0\": (0,0,0,255), "1": 0,0,255,255, "2": 0,255,255,255, "3": 255,255,255,255}"
    """
    logger.debug(
        "tile path %s, colormap: %s, stretch_range: %s, explicit_color_map: %s",
        keys,
        colormap,
        ast.literal_eval(stretch_range) if stretch_range else "",
        explicit_color_map,
    )
    try:
        source_db = _tile_db_from_keys(keys)
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
            if colormap == "explicit":
                # Check if we have an internal categorical colormap for this DB
                if source_db in CATEGORICAL_COLOR_MAPS.keys():
                    options["colormap"] = CATEGORICAL_COLOR_MAPS[source_db]
                elif not explicit_color_map:
                    raise MissingExplicitColourMapException()
                else:
                    # Use the provided categorical colormap
                    options["colormap"] = ast.literal_eval(explicit_color_map)
            else:
                options["colormap"] = colormap
        if stretch_range:
            options["stretch_range"] = ast.literal_eval(stretch_range)

        # Check the database exists
        if not _source_db_exists(source_db):
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
    except MissingExplicitColourMapException as err:
        handle_exception(logger, err)
        raise HTTPException(
            status_code=400,
            detail=f"colormap=explicit requires explicit_color_map to be included",
        )
    except SourceDBDoesNotExistException as err:
        handle_exception(logger, err)
        raise HTTPException(
            status_code=400,
            detail=f"source database {source_db} does not exist in tiles metastore",
        )
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)
