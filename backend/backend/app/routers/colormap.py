"""
Colormap
"""

import inspect
import json

from fastapi import APIRouter, HTTPException
from fastapi.logger import logger

from backend.app import schemas
from backend.app.internal.helpers import handle_exception


router = APIRouter(tags=["colormap"])


def _get_colormap(options: schemas.ColorMapOptions) -> schemas.ColorMap:
    """
    Retrieve colormap
    """
    from terracotta.handlers.colormap import colormap

    _colormap = colormap(**options.model_dump())
    return schemas.ColorMap.model_validate({"colormap": _colormap})


@router.get("", response_model=schemas.ColorMap)
async def get_colormap(
    colormap: str,
    stretch_range: str,
    num_values: int = 255,
) -> schemas.ColorMap:
    """
    Retrieve colormap values.  e.g. `colormap?colormap=reds&stretch_range=[0,10]`

    ::param colormap str The name of the colormap, e.g. `colormap=reds`

    ::param stretch_range iterable The url-encoded stretch-range over-which the colors should be generated

        e.g. `stretch_range=[0,10]`

    ::kwarg num_values int Number of values to generate in the colormap
    """
    try:
        stretch_range_arg: list[int] = json.loads(stretch_range)
        if len(stretch_range_arg) != 2:
            raise ValueError
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail="stretch_range should be a list of two values (min and max), like `stretch_range=[0,1]`",
        )

    logger.debug(
        "performing %s, using %s, %s, %s",
        inspect.stack()[0][3],
        colormap,
        stretch_range_arg,
        num_values,
    )
    try:

        options = schemas.ColorMapOptions(
            colormap=colormap,
            stretch_range=stretch_range_arg,
            num_values=num_values,
        )

        res = _get_colormap(options)

        return res
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)
