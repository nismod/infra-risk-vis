# Publishing the serving bundle (design)

Status: design only; deploy sync is not implemented.

irv-datapkg ADR 0007 leaves Web Mercator clipping, terracotta, the Zarr stack,
MBTiles and PostGIS in infra-risk-vis, and asks infra-risk-vis to publish its
serving bundle with notes on addressing it directly. This note describes that
bundle and how it would be published next to (not inside) the irv-datapkg
catalogue.

## What is in the bundle

| Product | Built by | Local path | Format |
| --- | --- | --- | --- |
| Raster tiles source | `cloud_optimise_raster` | `raster/cog/{dataset}/{key}.tif` | COG, EPSG:3857, clipped to `raster_bounds`, nearest resampling, overviews |
| Raster stack | `raster_dataset_tiffs_to_zarr` | `raster/stack.zarr/{dataset}/...` | Zarr v3, one group per dataset, dimensions from `metadata.json` `keys` |
| Vector tiles | tippecanoe (see README) | `tileserver/vector/data/*.mbtiles` | MBTiles (MVT) |
| Layer index | `layers.csv`, `metadata.json` | `pipelines/{dataset}/` | CSV / JSON |

PostGIS and the terracotta metadata databases are not published: they are
rebuilt from the bundle (`ingest_rasters`, `ingest_metadata`).

## Layout under a deploy prefix

One deploy prefix per deployment, e.g. `{DEPLOY} = s3://bucket/irv-serving/global`:

```
{DEPLOY}/
  {deploy_version}/                  # immutable, e.g. 2026.1-1 (datapkg release + deploy revision)
    manifest.json                    # what was built from what (below)
    datapkg.lock.yaml                # copy of etl/datapkg.lock.yaml used for the build
    raster/cog/{dataset}/{key}.tif
    raster/cog/{dataset}/layers.csv
    raster/cog/{dataset}/metadata.json
    raster/stack.zarr/
    vector/{name}.mbtiles
  current -> {deploy_version}        # a small JSON pointer object, not a symlink
```

- `deploy_version` directories are immutable, like irv-datapkg versions; the
  live site reads one fixed version, so a rollback is a pointer change.
- `manifest.json` lists every file with size and sha256, the infra-risk-vis
  commit, the irv-datapkg release and, per COG, the catalogue item it was
  derived from (`datapkg_item` from `layers.csv`). That gives provenance back to
  the upstream source through the STAC item.
- Licences and citations come from the STAC collections (copied into
  `metadata.json` by `scripts/metadata_from_catalogue.py`), so redistribution
  terms travel with each dataset.
- The deploy prefix is separate from the irv-datapkg `{BASE}`: these are
  visualisation products (Web Mercator, clipped at about 84 degrees), not T1
  data, and must not be mistaken for it.

A later sync step would be `aws s3 sync raster/cog {DEPLOY}/{deploy_version}/raster/cog`
(and the same for the Zarr store and MBTiles), then write `manifest.json` and
move `current`. The server's `tileserver/raster/data` volume can be filled from
the same prefix instead of from the build machine.

## Addressing it directly

GIS and web users can use the published files without the infra-risk-vis API:

- **COGs over HTTP.** Every COG supports range requests, so QGIS, GDAL and
  rasterio read windows without downloading the file:
  `gdalinfo /vsicurl/https://{host}/{deploy_version}/raster/cog/aqueduct/{key}.tif`,
  or in QGIS *Layer > Add Raster Layer > Protocol: HTTP(S)*. Browsers can
  render them with `geotiff.js` or OpenLayers' `GeoTIFF` source. They are
  EPSG:3857 and clipped, so for analysis prefer the T1 file in the irv-datapkg
  catalogue (source grid, not clipped), found from the `datapkg_item`.
- **Choosing a layer.** `layers.csv` beside the COGs maps key columns
  (hazard, rp, rcp, epoch, ...) to file names, as used by terracotta. The
  backend tile route keeps working for XYZ tiles:
  `/tiles/{domain}/{key values...}/{z}/{x}/{y}.png`, with domains and key order
  listed by `/tiles/sources` (paths relative to the backend API root).
- **Zarr stack.** `xarray.open_zarr("https://{host}/{deploy_version}/raster/stack.zarr", group="aqueduct")`
  gives the full n-dimensional array (lat/lon plus scenario dimensions) for
  notebooks and web clients such as `zarr.js`; chunks are fetched lazily.
- **MBTiles.** Download and open directly in QGIS, or serve with any MBTiles
  server (tileserver-gl, martin); the public site serves them at
  `/vector/data/{name}.json` (TileJSON) and `/vector/data/{name}/{z}/{x}/{y}.pbf`.
- **Pinning.** Address a fixed `{deploy_version}` for reproducible work; use
  `current` only for "whatever the site shows".

## Open questions

- Public HTTP host for the prefix (CDN in front of S3, CORS for browser reads).
- Whether to also publish a STAC catalogue for the bundle (items for each COG
  with `derived_from` links to irv-datapkg items) so STAC browsers can find it.
- Whether to publish PMTiles alongside MBTiles for direct browser use without
  a tile server.
