"""Scripted steps of the nbs-adaptation pipeline (pipelines/nbs-adaptation/Snakefile).

Replaces the extraction parts of the notebooks:

- extract-nbs-adaptation-boundaries.ipynb: boundaries and HydroBASINS to
  GeoJSONSeq polygons (with a ``bbox_wkt`` property) and centroid points, for
  tippecanoe;
- extract-nbs-adaptation-opportunities.ipynb, "Extract": the opportunity
  areas joined with EAD, one GeoParquet per NbS type, with every hazard column
  the later notebook cells use (missing ones zero-filled).

Inputs come from a published open-gira results run (layers
``admin_boundaries``, ``hydrobasins_gadm``, ``nbs_opportunities``) or, as a
fallback, from an open-gira results directory (``legacy_results_dir``).
"""

from pathlib import Path

ADMIN_COLUMNS = {
    0: ["GID_0", "NAME_0"],
    1: ["GID_0", "NAME_0", "GID_1", "NAME_1"],
    2: ["GID_0", "NAME_0", "GID_1", "NAME_1", "GID_2", "NAME_2"],
}

# raw_data/boundaries/{name}.geoparquet: (results layer, item filters,
# columns kept, MBTiles name)
BOUNDARIES = {
    **{
        f"admin-level-{n}": ("admin_boundaries", {"admin_level": n}, cols + ["geometry"], f"adm{n}")
        for n, cols in ADMIN_COLUMNS.items()
    },
    "hybas_lev12_v1c_with_gadm_codes": (
        "hydrobasins_gadm",
        {"level": 12},
        ["HYBAS_ID"] + ADMIN_COLUMNS[2] + ["geometry"],
        "hybas",
    ),
}

# where the boundaries are in an open-gira results directory
LEGACY_BOUNDARY_PATHS = {
    **{f"admin-level-{n}": f"input/admin-boundaries/admin-level-{n}.geoparquet" for n in ADMIN_COLUMNS},
    "hybas_lev12_v1c_with_gadm_codes": "input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet",
}

NBS_TYPES = ("river_basin_afforestation", "mangrove", "landslide_slope_vegetation")
AGGS = ("MAX", "MEAN", "MIN")
RCPS = ("rcp4p5", "rcp8p5")
EPOCHS = ("2030", "2050", "2080")
SECTORS = ("rail", "road")


def expected_columns(nbs: str) -> list:
    """Hazard columns the opportunities notebook reads for ``nbs`` (as listed
    in its "Extract" cell)."""
    cols = []
    if nbs == "river_basin_afforestation":
        cols += [f"hazard-floodMapGL_EAD__{s}" for s in SECTORS]
        cols += [f"hazard-inunriver_historical_{a}_1980_EAD__{s}" for a in AGGS for s in SECTORS]
        cols += [
            f"hazard-inunriver_{r}_{a}_{e}_EAD__{s}"
            for r in RCPS
            for a in AGGS
            for e in EPOCHS
            for s in SECTORS
        ]
    elif nbs == "mangrove":
        cols += [
            f"hazard-GFM_{dem}_{year}_EAD__{s}"
            for dem in ("MERITDEM1km", "NASADEM1km")
            for year in ("2018", "2050")
            for s in SECTORS
        ]
        cols += [
            f"hazard-inuncoast_historical_{a}_hist_{p}_EAD__{s}"
            for a in AGGS
            for p in ("0", "5")
            for s in SECTORS
        ]
        cols += [
            f"hazard-inuncoast_{r}_{a}_{e}_{p}_EAD__{s}"
            for r in RCPS
            for a in AGGS
            for e in EPOCHS
            for p in ("0", "5")
            for s in SECTORS
        ]
    elif nbs == "landslide_slope_vegetation":
        cols += [
            f"hazard-{h}__{s}_{v}_EAD__{s}"
            for s in SECTORS
            for h in ("ls_eq_tiled", "LS_RF_Median_1980-2018", "_landslide_sum")
            for v in ("damage_fraction", "lower", "upper")
        ]
    else:
        raise ValueError(f"unknown NbS type {nbs!r}: expected one of {NBS_TYPES}")
    return cols


def ensure_columns(gdf, columns):
    """Add any missing ``columns`` as zeros (as the notebook did per slice)."""
    for col in columns:
        if col not in gdf.columns:
            gdf[col] = 0
    return gdf


def _bbox_column(gdf):
    """In open-gira results runs, a GeoParquet 1.1 bbox covering column."""
    return gdf.drop(columns=["bbox"], errors="ignore")


def nbs_with_ead(sources, nbs: str, output) -> Path:
    """Combine opportunity files (a results-run item, or open-gira slices) into
    one GeoParquet with every expected hazard column."""
    import geopandas
    import pandas
    import pyarrow.parquet as pq

    dfs = []
    for path in sources:
        if "HYBAS_ID" not in pq.ParquetFile(path).schema_arrow.names:
            continue  # empty slice
        dfs.append(ensure_columns(_bbox_column(geopandas.read_parquet(path)), expected_columns(nbs)))
    if not dfs:
        raise ValueError(f"{nbs}: no non-empty opportunity files in {list(map(str, sources))}")
    gdf = geopandas.GeoDataFrame(pandas.concat(dfs), crs=dfs[0].crs)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(output)
    return output


def boundary_geojsonseq(src, name: str, polygons_out, points_out) -> None:
    """Boundaries to GeoJSONSeq polygons (plus ``bbox_wkt``) and centroid points."""
    import geopandas
    import shapely

    cols = BOUNDARIES[name][2]
    df = geopandas.read_parquet(src)
    if any(c not in df.columns for c in cols) and any(n in cols for n in df.index.names):
        df = df.reset_index()  # legacy HydroBASINS are indexed by HYBAS_ID
    df = df[cols].copy()
    bounds = df.bounds
    bbox = shapely.box(bounds.minx, bounds.miny, bounds.maxx, bounds.maxy)
    df["bbox_wkt"] = geopandas.GeoSeries(bbox, index=df.index).to_wkt()
    for out in (polygons_out, points_out):
        Path(out).parent.mkdir(parents=True, exist_ok=True)
    df.to_file(polygons_out, driver="GeoJSONSeq", engine="pyogrio")
    df.geometry = df.geometry.centroid
    df.to_file(points_out, driver="GeoJSONSeq", engine="pyogrio")
