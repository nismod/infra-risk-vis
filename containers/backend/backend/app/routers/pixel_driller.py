"""
Pixel Driller Service - Query Zarr raster data for pixel values at specific coordinates
"""

import logging
import traceback
from pathlib import Path
from typing import Optional, Union

import numpy as np
import xarray as xr
from fastapi import APIRouter, HTTPException
from fastapi.logger import logger
from pyproj import Transformer
from sqlalchemy import select

from backend.app import schemas
from backend.app.routers.tiles import _source_options, _tile_db_from_domain
from backend.config import ZARR_BASE_PATH
from backend.db import models
from backend.db.database import SessionDep

router = APIRouter(tags=["pixel-driller"])


def generate_layer_id(
    domain: str,
    keys: dict[str, Union[str, int, float]],
    key_order: list[str],
) -> str:
    """
    Generate a unique layer ID from domain and keys.
    
    Format: domain/value1/value2/... (matching Terracotta URL format)
    Uses Terracotta key order for consistency with tile URLs.
    """
    # Use Terracotta key order, fallback to sorted if key not found
    id_parts = [domain]
    for key in key_order:
        if key in keys:
            id_parts.append(str(keys[key]))
    return "/".join(id_parts)


def domain_to_zarr_group(domain: str, zarr_groups: dict) -> Optional[str]:
    """
    Map domain name (from database) to Zarr group prefix (pipeline directory name).
    
    The ETL uses pipeline directory names (e.g., 'gem_earthquake') as Zarr group prefixes,
    but the database uses the 'domain' from metadata.json (e.g., 'earthquake').
    
    Args:
        domain: Domain name from raster_tile_sources table
        zarr_groups: Dict of all Zarr group paths (from xr.open_groups())
    
    Returns:
        Zarr group prefix (dataset name) if found, None otherwise
    """
    # Known mappings where domain != pipeline directory name
    DOMAIN_TO_DATASET = {
        "earthquake": "gem_earthquake",
        "cyclone_storm": "storm",
        "cyclone_iris": "iris",
        "population": "ghsl_pop",
        "buildings": "ghsl_buildings",
        "landslide": "landslide_arup",
    }
    
    # First, try known mapping
    if domain in DOMAIN_TO_DATASET:
        dataset = DOMAIN_TO_DATASET[domain]
        # Check if this dataset exists in Zarr
        if any(g.startswith(f"/{dataset}/") for g in zarr_groups.keys()):
            return dataset
    
    # Second, try domain as-is (for cases where they match)
    if any(g.startswith(f"/{domain}/") for g in zarr_groups.keys()):
        return domain
    
    # Third, try to find a group that contains the domain name
    # (e.g., if domain is "earthquake", look for groups starting with "gem_earthquake")
    for group_path in zarr_groups.keys():
        # Extract the first part of the path (dataset name)
        if group_path.startswith("/"):
            parts = group_path[1:].split("/")
            if len(parts) > 0:
                dataset = parts[0]
                # Check if domain is a substring of dataset or vice versa
                if domain in dataset or dataset in domain:
                    return dataset
    
    return None


def query_zarr_pixel(
    zarr_store_path: Path,
    domain: str,
    lat: float,
    lon: float,
    session: SessionDep,
) -> list[schemas.PixelDrillerResult]:
    """
    Query Zarr store for pixel values at given coordinates.
    
    Uses Terracotta to discover valid key combinations, then queries Zarr using xarray.
    Coordinates are transformed from EPSG:4326 (degrees) to EPSG:3857 (Web Mercator meters)
    to match the Zarr store coordinate system.
    Similar approach to Jamaica implementation but adapted for multi-dimensional structure.
    """
    results = []
    
    logger.debug(f"[{domain}] Starting query_zarr_pixel for point ({lat}, {lon})")
    
    # Transform coordinates from EPSG:4326 (WGS84 degrees) to EPSG:3857 (Web Mercator meters)
    # Zarr store uses EPSG:3857 but coordinates are named 'lat'/'lon' (misleading naming)
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    lon_proj, lat_proj = transformer.transform(lon, lat)
    logger.debug(f"[{domain}] Transformed coordinates: ({lat}, {lon}) -> ({lat_proj}, {lon_proj})")
    
    try:
        # OPTIMIZATION: Check Zarr groups first before querying Terracotta (which is slow)
        # This avoids expensive Terracotta queries for domains that don't have Zarr data
        zarr_path_str = str(zarr_store_path)
        if not zarr_path_str or zarr_path_str.strip() == '':
            logger.warning(f"[{domain}] Empty Zarr store path. Check ZARR_BASE_PATH environment variable.")
            return results
        
        logger.debug(f"[{domain}] Opening Zarr store to check for domain groups...")
        try:
            # xarray.open_groups() returns a dict-like object with all group paths as keys
            # Groups are like '/aqueduct/depth_coastal', '/dem/elevation', etc.
            zarr_groups = xr.open_groups(zarr_path_str)
            logger.debug(f"[{domain}] Successfully opened Zarr store, found {len(zarr_groups)} total groups")
        except Exception as e:
            logger.warning(f"[{domain}] Failed to open Zarr store at '{zarr_path_str}': {e}")
            return results
        
        # Map domain to Zarr group prefix (dataset name)
        zarr_dataset = domain_to_zarr_group(domain, zarr_groups)
        if not zarr_dataset:
            logger.debug(f"[{domain}] No matching Zarr group found for domain {domain}, skipping Terracotta query")
            return results
        
        logger.debug(f"[{domain}] Mapped domain '{domain}' to Zarr dataset '{zarr_dataset}', proceeding with Terracotta query")
        
        # Now that we know the domain exists in Zarr, fetch valid combinations from Terracotta
        # Get Terracotta database name for this domain
        logger.debug(f"[{domain}] Getting Terracotta database name")
        source_db = _tile_db_from_domain(domain)
        logger.debug(f"[{domain}] Terracotta database: {source_db}")
        
        # Get all valid key combinations from Terracotta
        # Note: _source_options creates a driver internally, but we avoid creating a second one
        logger.debug(f"[{domain}] Fetching valid combinations from Terracotta...")
        valid_combinations = _source_options(source_db)
        logger.debug(f"[{domain}] Found {len(valid_combinations)} valid combinations")
        
        if not valid_combinations:
            logger.debug(f"[{domain}] No valid combinations found for domain {domain}")
            return results
        
        # Extract key order from the first combination
        # In Python 3.7+, dicts maintain insertion order, and Terracotta returns keys in order
        # This avoids creating a second driver connection just to get key order
        key_order = list(valid_combinations[0].keys())
        logger.debug(f"[{domain}] Key order: {key_order}")
        
        # Find all groups that match {zarr_dataset}/{var} pattern
        var_names = []
        dataset_prefix = f"/{zarr_dataset}/"
        logger.debug(f"[{domain}] Searching for groups matching prefix: {dataset_prefix}")
        for group_path in zarr_groups.keys():
            # Group paths from xarray.open_groups() start with '/'
            if group_path.startswith(dataset_prefix):
                # Extract var name (everything after dataset/)
                var_name = group_path[len(dataset_prefix):]
                # Only take direct children (not nested groups like dataset/var/subgroup)
                if '/' not in var_name:
                    var_names.append(var_name)
        
        logger.debug(f"[{domain}] Found {len(var_names)} var groups: {var_names}")
        if not var_names:
            logger.debug(f"[{domain}] No var groups found for domain {domain} (dataset {zarr_dataset}) in Zarr store")
            return results
        
        # Process each var group
        for var_idx, var in enumerate(var_names, 1):
            var_group_path = f"{zarr_dataset}/{var}"
            logger.debug(f"[{domain}] Processing var {var_idx}/{len(var_names)}: {var_group_path}")
            
            try:
                # Open the var group as xarray dataset
                logger.debug(f"[{domain}/{var}] Opening xarray dataset...")
                var_ds = xr.open_zarr(str(zarr_store_path), group=var_group_path)
                logger.debug(f"[{domain}/{var}] Dataset opened successfully")
            except (KeyError, OSError) as e:
                logger.warning(f"[{domain}/{var}] Failed to open {var_group_path}: {e}")
                continue
            
            # Get the data variable (named after the var)
            if var not in var_ds.data_vars:
                logger.warning(f"[{domain}/{var}] Data variable {var} not found in {var_group_path}")
                continue
            
            data_var = var_ds[var]
            logger.debug(f"[{domain}/{var}] Data variable shape: {data_var.shape}, dims: {data_var.dims}")
            
            # Get dimension names
            dim_names = list(data_var.dims)
            non_spatial_dims = [d for d in dim_names if d not in ['lat', 'lon']]
            logger.debug(f"[{domain}/{var}] Non-spatial dimensions: {non_spatial_dims}")
            
            # Check bounds before querying
            if 'lat' not in var_ds.coords or 'lon' not in var_ds.coords:
                logger.warning(f"[{domain}/{var}] Coordinate arrays not found in {var_group_path}")
                continue
            
            logger.debug(f"[{domain}/{var}] Computing coordinate bounds...")
            lat_coords = var_ds.coords['lat']
            lon_coords = var_ds.coords['lon']
            lat_min, lat_max = float(lat_coords.min()), float(lat_coords.max())
            lon_min, lon_max = float(lon_coords.min()), float(lon_coords.max())
            logger.debug(f"[{domain}/{var}] Bounds: lat=[{lat_min}, {lat_max}], lon=[{lon_min}, {lon_max}]")
            
            # Check if point is within bounds (using projected coordinates)
            if lat_proj < lat_min or lat_proj > lat_max or lon_proj < lon_min or lon_proj > lon_max:
                logger.debug(f"[{domain}/{var}] Point ({lat_proj}, {lon_proj}) outside bounds, skipping")
                continue
            
            logger.debug(f"[{domain}/{var}] Point is within bounds, processing {len(valid_combinations)} combinations...")
            
            # OPTIMIZATION 1: Early filtering - pre-compute valid dimension values to filter out invalid combinations
            # This avoids expensive xarray operations on combinations that will fail anyway
            valid_dim_values = {}
            dim_dtypes = {}  # Store dtype info for type conversion during filtering
            if non_spatial_dims:
                logger.debug(f"[{domain}/{var}] Pre-computing valid dimension values for early filtering...")
                for dim in non_spatial_dims:
                    try:
                        dim_coords = var_ds.coords[dim]
                        dim_dtypes[dim] = dim_coords.dtype.kind  # Store dtype for later use
                        # Convert to numpy array (handles lazy arrays) and create a set for fast lookup
                        # Use .to_numpy() or .values depending on xarray version
                        try:
                            coord_values = dim_coords.to_numpy()
                        except AttributeError:
                            # Fallback for older xarray versions
                            coord_values = dim_coords.values
                        
                        # Convert to appropriate type and create a set for fast lookup
                        if dim_coords.dtype.kind == 'U':  # String array
                            valid_dim_values[dim] = set(str(v) for v in coord_values)
                        elif dim_coords.dtype.kind in ['i', 'u']:  # Integer array
                            valid_dim_values[dim] = set(int(v) for v in coord_values)
                        else:  # Float array
                            valid_dim_values[dim] = set(float(v) for v in coord_values)
                        logger.debug(f"[{domain}/{var}] Dimension {dim} has {len(valid_dim_values[dim])} valid values")
                    except Exception as e:
                        logger.warning(f"[{domain}/{var}] Failed to pre-compute values for dimension {dim}: {e}")
                        # If we can't pre-compute, skip early filtering for this dimension
                        # The combination will be checked during actual selection
                        valid_dim_values[dim] = None
            
            # Filter combinations early - only keep those that match var dimensions AND have valid values
            filtered_combinations = []
            for combo in valid_combinations:
                # Check if keys match
                combo_keys = set(combo.keys())
                var_dims = set(non_spatial_dims)
                if combo_keys != var_dims:
                    continue
                
                # Check if all dimension values exist in Zarr coordinates
                is_valid = True
                for dim in non_spatial_dims:
                    target_value = combo[dim]
                    # Convert to appropriate type for comparison
                    if dim in valid_dim_values and valid_dim_values[dim] is not None:
                        dtype_kind = dim_dtypes[dim]
                        if dtype_kind == 'U':  # String
                            target_value = str(target_value)
                        elif dtype_kind in ['i', 'u']:  # Integer
                            try:
                                target_value = int(target_value)
                            except (ValueError, TypeError):
                                target_value = str(target_value)
                        else:  # Float
                            try:
                                target_value = float(target_value)
                            except (ValueError, TypeError):
                                target_value = str(target_value)
                        
                        if target_value not in valid_dim_values[dim]:
                            is_valid = False
                            break
                    # If valid_dim_values[dim] is None, we couldn't pre-compute, so skip early filtering
                    # The combination will be checked during actual selection
                
                if is_valid:
                    filtered_combinations.append(combo)
            
            logger.debug(f"[{domain}/{var}] Filtered from {len(valid_combinations)} to {len(filtered_combinations)} valid combinations")
            
            if not filtered_combinations:
                logger.debug(f"[{domain}/{var}] No valid combinations after filtering, skipping")
                continue
            
            # OPTIMIZATION 2: Select spatial coordinates first (once per var, not per combination)
            # This gives us a spatial slice that we can then filter by dimensions
            # This is much more efficient since spatial selection happens once instead of N times
            spatial_sel = {"lat": lat_proj, "lon": lon_proj}
            try:
                # Select spatial coordinates first - this is the same for all combinations
                spatial_slice = data_var.sel(spatial_sel, method="nearest")
                logger.debug(f"[{domain}/{var}] Spatial slice shape: {spatial_slice.shape}, dims: {spatial_slice.dims}")
            except Exception as e:
                logger.warning(f"[{domain}/{var}] Failed to select spatial coordinates: {e}")
                continue
            
            # Process each valid combination from Terracotta
            combo_count = 0
            for combo_idx, combo in enumerate(filtered_combinations, 1):
                if combo_idx % 100 == 0:
                    logger.debug(f"[{domain}/{var}] Processed {combo_idx}/{len(filtered_combinations)} combinations...")
                
                combo_count += 1
                
                # Build selection dictionary for non-spatial dimensions only
                non_spatial_sel = {}
                for dim in non_spatial_dims:
                    target_value = combo[dim]
                    # Convert to appropriate type for matching
                    dim_coords = var_ds.coords[dim]
                    if dim_coords.dtype.kind == 'U':  # String array
                        target_value = str(target_value)
                    elif dim_coords.dtype.kind in ['i', 'u']:  # Integer array
                        try:
                            target_value = int(target_value)
                        except (ValueError, TypeError):
                            target_value = str(target_value)
                    else:  # Float array
                        try:
                            target_value = float(target_value)
                        except (ValueError, TypeError):
                            target_value = str(target_value)
                    
                    non_spatial_sel[dim] = target_value
                
                try:
                    # Now select non-spatial dimensions on the already-filtered spatial slice
                    # This is much faster since we're working with a smaller array
                    if non_spatial_sel:
                        selected = spatial_slice.sel(non_spatial_sel)
                    else:
                        selected = spatial_slice
                    
                    # logger.debug(f"[{domain}/{var}] Selection successful, shape: {selected.shape}")
                    
                    # Extract scalar value
                    value = float(selected.values.item()) if selected.size == 1 else None
                    
                    # Handle NaN/None
                    if value is not None and (np.isnan(value) or np.isinf(value)):
                        value = None
                    
                    # Generate layer ID (using Terracotta key order for consistency)
                    layer_id = generate_layer_id(domain, combo, key_order)
                    
                    results.append(
                        schemas.PixelDrillerResult(
                            value=value,
                            layer=schemas.PixelDrillerLayer(
                                domain=domain,
                                keys=combo,
                                id=layer_id,
                            ),
                        )
                    )
                    # logger.debug(f"[{domain}/{var}] Added result for combo {combo_idx}, total results: {len(results)}")
                except (KeyError, ValueError) as e:
                    # This should rarely happen now due to early filtering, but keep for safety
                    logger.debug(f"[{domain}/{var}] Selection failed for combo {combo_idx} ({combo}): {e}")
                    continue
                except Exception as e:
                    logger.warning(f"[{domain}/{var}] Unexpected error processing combo {combo_idx} ({combo}): {e}")
                    continue
            
            logger.debug(f"[{domain}/{var}] Processed {combo_count} matching combinations, added {len([r for r in results if r.layer.domain == domain])} results")
        
    except KeyError as e:
        logger.warning(f"[{domain}] Domain or component not found: {e}")
        logger.debug(f"[{domain}] KeyError traceback:\n{traceback.format_exc()}")
    except Exception as e:
        # Log error but don't fail the entire request - skip this domain
        # Reduce log noise for connection errors (they're expected when DB is exhausted)
        if "connection" in str(e).lower() or "database" in str(e).lower():
            logger.warning(f"[{domain}] Database connection error: {e}")
        else:
            logger.error(f"[{domain}] Error querying Zarr: {e}")
            logger.debug(f"[{domain}] Full traceback:\n{traceback.format_exc()}")
    
    logger.debug(f"[{domain}] Completed query_zarr_pixel, returning {len(results)} results")
    return results


@router.get("/point/{lon}/{lat}", response_model=schemas.PixelDrillerResponse)
def get_pixel_values(
    lon: float,
    lat: float,
    session: SessionDep,
) -> schemas.PixelDrillerResponse:
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
        raise HTTPException(status_code=400, detail="Longitude must be between -180 and 180")
    if not (-90 <= lat <= 90):
        raise HTTPException(status_code=400, detail="Latitude must be between -90 and 90")
    
    # Validate Zarr store path
    if not ZARR_BASE_PATH or not ZARR_BASE_PATH.strip():
        raise HTTPException(
            status_code=500,
            detail="ZARR_BASE_PATH environment variable is not set or is empty",
        )
    
    zarr_path = Path(ZARR_BASE_PATH)
    if not zarr_path.exists():
        raise HTTPException(
            status_code=500,
            detail=f"Zarr store not found at {ZARR_BASE_PATH}",
        )
    
    logger.debug(f"Using Zarr store at: {zarr_path}")
    
    # Get all domains from RasterTileSource
    logger.debug("Querying RasterTileSource for domains...")
    try:
        all_domains = session.scalars(select(models.RasterTileSource)).all()
        domain_list = [d.domain for d in all_domains]
        logger.debug(f"Found {len(domain_list)} domains: {domain_list}")
    except Exception as e:
        logger.error(f"Failed to query RasterTileSource: {e}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to query domains: {str(e)}",
        )
    
    # TEMPORARY: Filter to only test domains for faster testing
    TEST_DOMAINS = {"aqueduct", "cyclone_storm", "landslide", "isimip"}
    domains = [d for d in all_domains if d.domain in TEST_DOMAINS]
    logger.debug(f"Filtered to {len(domains)} test domains: {[d.domain for d in domains]}")
    
    # Query Zarr for each domain
    # Wrap in try/except to handle connection errors gracefully
    all_results = []
    for domain_idx, domain_meta in enumerate(domains, 1):
        logger.debug(f"Processing domain {domain_idx}/{len(domains)}: {domain_meta.domain}")
        try:
            domain_results = query_zarr_pixel(
                zarr_path,
                domain_meta.domain,
                lat,
                lon,
                session,
            )
            logger.debug(f"Domain {domain_meta.domain} returned {len(domain_results)} results")
            all_results.extend(domain_results)
        except Exception as e:
            # Skip domains that fail (e.g., connection errors, missing Zarr data)
            # Log with full traceback for debugging
            logger.error(f"Skipping domain {domain_meta.domain} due to error: {e}")
            logger.debug(f"Full traceback for domain {domain_meta.domain}:\n{traceback.format_exc()}")
            continue
    
    logger.debug(f"Completed processing all domains, total results: {len(all_results)}")
    
    try:
        return schemas.PixelDrillerResponse(
            point={"lat": lat, "lon": lon},
            results=all_results,
        )
    except Exception as e:
        logger.error(f"Failed to create response: {e}")
        logger.debug(f"Response creation traceback:\n{traceback.format_exc()}")
        raise HTTPException(
            status_code=500,
            detail=f"Failed to create response: {str(e)}",
        )

