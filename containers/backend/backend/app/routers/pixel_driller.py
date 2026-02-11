"""
Pixel Driller Service - Query Zarr raster data for pixel values at specific coordinates
"""

import itertools
import traceback

import numpy as np
import xarray as xr
from fastapi import APIRouter, HTTPException
from fastapi.logger import logger
from fastapi.responses import ORJSONResponse
from pyproj import Transformer

from backend.app import schemas
from backend.config import ZARR_GROUPS, ZARR_STORE

router = APIRouter(tags=["pixel-driller"])


def zarr_group_to_domain(group: str) -> str:
    """
    Map Zarr group prefix (pipeline directory name) to domain key (from
    terracotta database).

    The ETL uses pipeline directory names (e.g., 'gem_earthquake') as Zarr group
    prefixes, but the database uses the 'domain' from metadata.json (e.g.,
    'earthquake') in the raster_tile_sources table.

    Args:
        group: Zarr group prefix (dataset name)

    Returns:
        Domain key
    """
    group = group.split("/")[1]
    # Known mappings where domain != pipeline directory name
    GROUP_TO_DOMAIN = {
        "gem_earthquake": "earthquake",
        "storm": "cyclone_storm",
        "iris": "cyclone_iris",
        "ghsl_pop": "population",
        "ghsl_buildings": "buildings",
        "landslide_arup": "landslide",
    }
    try:
        domain = GROUP_TO_DOMAIN[group]
    except KeyError:
        domain = group

    return domain


def check_bounds(ds, lat_proj, lon_proj):
    # Check bounds before querying
    if "lat" not in ds.coords or "lon" not in ds.coords:
        return False

    lat_coords = ds.coords["lat"]
    lon_coords = ds.coords["lon"]
    lat_min, lat_max = float(lat_coords.min()), float(lat_coords.max())
    lon_min, lon_max = float(lon_coords.min()), float(lon_coords.max())

    # Check if point is within bounds (using projected coordinates)
    return not (
        lat_proj < lat_min
        or lat_proj > lat_max
        or lon_proj < lon_min
        or lon_proj > lon_max
    )


def query_zarr_pixel(
    ds: xr.Dataset,
    group: str,
    lat: float,
    lon: float,
) -> list[schemas.PixelDrillerResult]:
    """
    Query Zarr store for pixel values at given coordinates.

    Queries Zarr using xarray, relies on Zarr metadata for response.
    Coordinates are transformed from EPSG:4326 (degrees) to EPSG:3857 (Web Mercator meters)
    to match the Zarr store coordinate system.
    """
    logger.debug(f"[{group}] Starting query_zarr_pixel for point ({lat}, {lon})")

    # Transform coordinates from EPSG:4326 (WGS84 degrees) to EPSG:3857 (Web Mercator meters)
    # Zarr store uses EPSG:3857 but coordinates are named 'lat'/'lon' (misleading naming)
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    lon_proj, lat_proj = transformer.transform(lon, lat)
    if not check_bounds(ds, lat_proj, lon_proj):
        logger.debug(
            f"[{group}] Point({lat}, {lon}) -> Point({lat_proj}, {lon_proj}) outside bounds, skipping"
        )
        return []

    variable_names = list(ds.keys())
    spatial_sel = {"lat": lat_proj, "lon": lon_proj}
    try:
        # Select spatial coordinates first - this is the same for all combinations
        spatial_slice = ds.sel(spatial_sel, method="nearest").compute()
    except Exception as e:
        logger.warning(f"[{group}] Failed to select spatial coordinates: {e}")
        return []

    results = []
    for var in variable_names:
        df = spatial_slice[var].to_series().reset_index()
        for dict_ in df.to_dict(orient="records"):
            value = dict_[var]

            # Handle NaN/None
            if value is not None and (np.isnan(value) or np.isinf(value)):
                value = None
            del dict_[var]
            results.append(
                {
                    "value": value,
                    "layer": {
                        "domain": zarr_group_to_domain(group),
                        "keys": dict_,
                    },
                }
            )
    logger.debug(
        f"[{group}] Completed query_zarr_pixel, returning {len(results)} results"
    )
    return results


@router.get(
    "/point/{lon}/{lat}",
    # response_model=schemas.PixelDrillerResponse,
    response_class=ORJSONResponse,
)
def get_pixel_values(
    lon: float,
    lat: float,
):
    """
    Get pixel values for all layers at a specific lat/lon coordinate.

    Args:
        lon: Longitude in degrees (-180 to 180)
        lat: Latitude in degrees (-90 to 90)

    Returns:
        PixelDrillerResponse with point coordinates, count, and results array
    """
    # Validate coordinates
    if not (-180 <= lon <= 180):
        raise HTTPException(
            status_code=400, detail="Longitude must be between -180 and 180"
        )
    if not (-90 <= lat <= 90):
        raise HTTPException(
            status_code=400, detail="Latitude must be between -90 and 90"
        )

    groups = ZARR_GROUPS
    store = ZARR_STORE

    # Query Zarr for each domain
    # Wrap in try/except to handle connection errors gracefully
    all_results = []
    for group_idx, group in enumerate(groups, 1):
        logger.info(f"Processing group {group_idx}/{len(groups)}: {group}")
        try:
            group_results = query_zarr_pixel(
                store[group],
                group,
                lat,
                lon,
            )
            logger.debug(f"Domain {group} returned {len(group_results)} results")
            all_results.append(group_results)
        except Exception as e:
            # Skip groups that fail (e.g., connection errors, missing Zarr data)
            # Log with full traceback for debugging
            logger.error(f"Skipping group {group} due to error: {e}")
            logger.debug(f"Full traceback for group {group}:\n{traceback.format_exc()}")
            continue

    logger.debug(f"Completed processing all groups, total results: {len(all_results)}")

    try:
        return ORJSONResponse(
            {
                "point": {"lat": lat, "lon": lon},
                "results": list(itertools.chain.from_iterable(all_results)),
            }
        )
    except Exception as e:
        logger.error(f"Failed to create response: {e}")
        logger.debug(f"Response creation traceback:\n{traceback.format_exc()}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to create response: {str(e)}",
        )
