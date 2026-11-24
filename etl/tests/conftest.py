"""Fixtures: a tiny synthetic irv-datapkg catalogue built with the producer API.

Nothing is downloaded: every raster is a 20x10 COG written locally. Item
properties mimic what WP5 is expected to publish for aqueduct_flood and
storm_fixed_rp, so the real pipelines/*/datapkg.yml mappings are exercised.
"""

import shutil
import sys
from pathlib import Path

import pandas as pd
import pytest

ETL = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ETL / "scripts"))

catalog = pytest.importorskip("irv_datapkg.catalog")
np = pytest.importorskip("numpy")
rasterio = pytest.importorskip("rasterio")

COMMIT = "0123456789abcdef0123456789abcdef01234567"

# (item id, properties, layers.csv filename or None if not in layers.csv)
AQUEDUCT_ITEMS = [
    (
        "inuncoast_historical_wtsub_hist_rp0001_5",
        {"hazard": "coastal", "rcp": "baseline", "epoch": "baseline", "gcm": None, "return_period": 1},
        "inuncoast_historical_wtsub_hist_rp0001_5.tif",
    ),
    (
        "inunriver_rcp4p5_noresm1-m_2030_rp00100",
        {"hazard": "river", "rcp": "4.5", "epoch": 2030, "gcm": "NorESM1-M", "return_period": 100},
        "inunriver_rcp4p5_00000NorESM1-M_2030_rp00100.tif",
    ),
    (
        "inunriver_historical_000000000watch_1980_rp00010",
        {"hazard": "river", "rcp": "baseline", "epoch": "baseline", "gcm": "WATCH", "return_period": 10},
        "inunriver_historical_000000000WATCH_1980_rp00010.tif",
    ),
    (
        "inunriver_rcp8p5_noresm1-m_2080_rp01000",
        {"hazard": "river", "rcp": "8.5", "epoch": 2080, "gcm": "NorESM1-M", "return_period": 1000},
        None,  # published but not in the (subset) layers.csv
    ),
]

STORM_ITEMS = [
    ("present_rp10", {"return_period": 10, "gcm": "constant", "rcp": "baseline", "epoch": "baseline"}),
    ("cmcc_2050_rp10", {"return_period": 10, "gcm": "CMCC-CM2-VHR4", "rcp": "8.5", "epoch": 2050}),
]


def make_cog(path, value):
    from rasterio.io import MemoryFile
    from rasterio.transform import from_origin
    import rasterio.shutil

    data = np.full((1, 10, 20), value, dtype="float32")
    profile = dict(driver="GTiff", width=20, height=10, count=1, dtype="float32",
                   crs="EPSG:4326", transform=from_origin(0, 1, 0.1, 0.1), nodata=-9999.0)
    with MemoryFile() as mem:
        with mem.open(**profile) as ds:
            ds.write(data)
        with mem.open() as ds:
            rasterio.shutil.copy(ds, str(path), driver="COG", compress="DEFLATE")
    return path


def _collection(dataset, dims, title):
    return catalog.build_collection(
        dataset,
        title=title,
        description=f"Synthetic {title}",
        license="CC-BY-4.0",
        citation=f"Nobody (2026) Synthetic {title}.",
        doi="10.0000/synthetic." + dataset,
        upstream_url="https://example.org/" + dataset,
        upstream_version="v1",
        processing=["convert to COG"],
        dimensions=dims,
    )


@pytest.fixture(scope="session")
def catalogue(tmp_path_factory):
    """Synthetic catalogue; returns its base path. Release 2026.1 pins both datasets."""
    root = tmp_path_factory.mktemp("datapkg")
    src = root / "src"
    src.mkdir()
    base = root / "base"
    writer = catalog.CatalogWriter(base)

    writer.add_collection(_collection(
        "aqueduct_flood", ["hazard", "rcp", "epoch", "gcm", "return_period"], "aqueduct"))
    for n, (item_id, props, _) in enumerate(AQUEDUCT_ITEMS):
        path = make_cog(src / f"{item_id}.tif", float(n + 1))
        item = catalog.build_item(path, "aqueduct_flood", "v2-r1", item_id, properties=props,
                                  source_url="https://example.org/a.tif", retrieved="2026-01-15",
                                  commit=COMMIT)
        writer.add_item(item, path)

    writer.add_collection(_collection(
        "storm_fixed_rp", ["return_period", "gcm", "rcp", "epoch"], "storm"))
    for n, (item_id, props) in enumerate(STORM_ITEMS):
        path = make_cog(src / f"storm_{item_id}.tif", 10.0 + n)
        item = catalog.build_item(path, "storm_fixed_rp", "1.0.0-r1", item_id, properties=props,
                                  source_url="https://example.org/s.zip", retrieved="2026-01-15",
                                  commit=COMMIT)
        writer.add_item(item, path)

    writer.release("2026.1")
    return base


@pytest.fixture
def etl_root(tmp_path):
    """A scratch copy of etl/pipelines with a 3-row aqueduct layers.csv subset."""
    root = tmp_path / "etl"
    shutil.copytree(ETL / "pipelines", root / "pipelines",
                    ignore=shutil.ignore_patterns("*.ipynb", "*.png", "standalone"))
    layers = pd.read_csv(ETL / "pipelines/aqueduct/layers.csv", keep_default_na=False)
    wanted = [f for _, _, f in AQUEDUCT_ITEMS if f]
    subset = layers[layers.filename.isin(wanted)]
    assert len(subset) == len(wanted)
    subset.to_csv(root / "pipelines/aqueduct/layers.csv", index=False)
    storm = pd.read_csv(ETL / "pipelines/storm/layers.csv", keep_default_na=False)
    storm = storm[storm.filename.isin([
        "STORM_FIXED_RETURN_PERIODS_constant_10_YR_RP.tif",
        "STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_10_YR_RP.tif",
    ])]
    assert len(storm) == 2
    storm.to_csv(root / "pipelines/storm/layers.csv", index=False)
    return root


@pytest.fixture
def cache(tmp_path):
    return str(tmp_path / "cache")
