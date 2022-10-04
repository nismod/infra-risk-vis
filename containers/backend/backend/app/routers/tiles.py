"""
Tile Service wrapping Terracotta Python API
"""
from sys import getsizeof
from typing import BinaryIO, List, Tuple, Union
import ast

from fastapi import APIRouter, Depends
from fastapi.logger import logger
from fastapi.responses import FileResponse
from starlette.responses import StreamingResponse
from mercantile import tile
from sqlalchemy import desc
from sqlalchemy.orm import Session
from geoalchemy2 import functions


from app import schemas
from app.dependencies import get_db
from db import models


router = APIRouter(tags=["tiles"])


def _get_singleband_image(
    keys: str, tile_xyz: Tuple[int, int, int] = None, options: dict = {}
) -> BinaryIO:
    from app.internal.tiles.singleband import singleband

    parsed_keys = [key for key in keys.split("/") if key]

    if options.get("colormap", "") == "explicit":
        options["colormap"] = options.pop("explicit_color_map")
    logger.debug(
        "parsed_keys: %s, tile_xyz: %s, options: %s", parsed_keys, tile_xyz, options
    )
    return singleband(parsed_keys, tile_xyz=tile_xyz, **options)


@router.get(
    "/{keys:path}/{tile_z:int}/{tile_x:int}/{tile_y:int}.png",
    response_class=FileResponse,
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

    # Check the keys are appropriate for the given type
    options = {}
    if colormap:
        options["colormap"] = colormap
    if stretch_range:
        options["stretch_range"] = ast.literal_eval(stretch_range)

    # Generate the tile
    image = _get_singleband_image(keys, [tile_x, tile_y, tile_z], options)
    logger.debug("tile image of size returned: %s, %s", getsizeof(image), type(image))

    # Return the tile as stream
    return StreamingResponse(image, media_type="image/png")
