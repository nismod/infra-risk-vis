"""
Colormap
"""
import ast
from typing import List, Optional, Union
import inspect

from fastapi import APIRouter, HTTPException
from fastapi.logger import logger


from app import schemas
from app.internal.helpers import handle_exception
from app.exceptions import SourceDBDoesNotExistException, DomainAlreadyExistsException


router = APIRouter(tags=["colormap"])


def _get_colormap(options: schemas.ColorMapOptions) -> schemas.ColorMap:
    """
    Retrieve colormap
    """
    from terracotta.handlers.colormap import colormap

    _colormap = colormap(**options.dict())
    return schemas.ColorMap.parse_obj({"colormap": _colormap})


@router.get("", response_model=schemas.ColorMap)
async def get_colormap(
    colormap: str,
    stretch_range: Union[str, None],
    num_values: int = 255,
) -> schemas.ColorMap:
    """
    Retrieve colormap
    """
    logger.debug(
        "performing %s, using %s, %s, %s",
        inspect.stack()[0][3],
        colormap,
        ast.literal_eval(stretch_range),
        num_values,
    )
    try:

        options = schemas.ColorMapOptions(
            colormap=colormap,
            stretch_range=ast.literal_eval(stretch_range),
            num_values=num_values,
        )

        res = _get_colormap(options)

        return res
    except Exception as err:
        handle_exception(logger, err)
        raise HTTPException(status_code=500)
