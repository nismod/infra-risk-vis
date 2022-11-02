"""handlers/singleband.py
Handle /singleband API endpoint.

Taken directly from:  https://github.com/DHI-GRAS/terracotta/blob/main/LICENSE

MIT License

Copyright (c) 2018 DHI GRAS

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from typing import (
    List,
    OrderedDict,
    Sequence,
    Mapping,
    Union,
    Tuple,
    Optional,
    TypeVar,
    cast,
    BinaryIO,
)

import collections

from terracotta import get_settings, get_driver, image, xyz
from terracotta.profile import trace

Number = TypeVar("Number", int, float)
RGBA = Tuple[Number, Number, Number, Number]


def database_keys(tc_driver_path: str) -> OrderedDict:
    """Keys and their orderinng in the given DB"""
    settings = get_settings()
    driver = get_driver(tc_driver_path, provider=settings.DRIVER_PROVIDER)
    return driver.get_keys()


def all_datasets(tc_driver_path: str) -> dict:
    """All datasets in the given database"""
    settings = get_settings()
    driver = get_driver(tc_driver_path, provider=settings.DRIVER_PROVIDER)
    return driver.get_datasets()


def singleband(
    tc_driver_path: str,
    keys: Union[Sequence[str], Mapping[str, str]],
    tile_xyz: Tuple[int, int, int] = None,
    *,
    colormap: Union[str, Mapping[Number, RGBA], None] = None,
    stretch_range: Tuple[Number, Number] = None,
    tile_size: Tuple[int, int] = None
) -> BinaryIO:
    """Return singleband image as PNG"""

    cmap_or_palette: Union[str, Sequence[RGBA], None]

    if stretch_range is None:
        stretch_min, stretch_max = None, None
    else:
        stretch_min, stretch_max = stretch_range

    preserve_values = isinstance(colormap, collections.abc.Mapping)

    settings = get_settings()
    if tile_size is None:
        tile_size = settings.DEFAULT_TILE_SIZE

    driver = get_driver(tc_driver_path, provider=settings.DRIVER_PROVIDER)

    with driver.connect():
        metadata = driver.get_metadata(keys)
        tile_data = xyz.get_tile_data(
            driver, keys, tile_xyz, tile_size=tile_size, preserve_values=preserve_values
        )

    if preserve_values:
        # bin output image into supplied labels, starting at 1
        colormap = cast(Mapping, colormap)

        labels, label_colors = list(colormap.keys()), list(colormap.values())

        cmap_or_palette = label_colors
        out = image.label(tile_data, labels)
    else:
        # determine stretch range from metadata and arguments
        stretch_range_ = list(metadata["range"])

        if stretch_min is not None:
            stretch_range_[0] = stretch_min

        if stretch_max is not None:
            stretch_range_[1] = stretch_max

        cmap_or_palette = cast(Optional[str], colormap)
        out = image.to_uint8(tile_data, *stretch_range_)

    return image.array_to_png(out, colormap=cmap_or_palette)
