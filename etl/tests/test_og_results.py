"""Vector layers from a published open-gira results run (scripts/og_results.py,
scripts/nbs_adaptation.py and the osm_roads, osm_rail and nbs-adaptation
Snakefiles). Everything runs against a synthetic run; nothing is downloaded.
"""

import json
import os
import re
from collections import Counter
from pathlib import Path

import pytest

from conftest import ETL

gpd = pytest.importorskip("geopandas")
pytest.importorskip("irv_datapkg.catalog")

import nbs_adaptation as nbs  # noqa: E402
import og_results  # noqa: E402
import og_results_fixture as fx  # noqa: E402
from test_snakemake import jobs, snakemake  # noqa: E402


@pytest.fixture(scope="module")
def base(tmp_path_factory):
    return fx.build(tmp_path_factory.mktemp("og"))


def cfg(base, **kw):
    return og_results.results_config({"open_gira_results": {"base": str(base), "run": fx.RUN, **kw}})


def test_config_and_activation(base, tmp_path):
    assert not og_results.active(og_results.results_config({}))
    assert not og_results.active(og_results.results_config(
        {"open_gira_results": {"base": "s3://x", "run": "None"}}))  # --config null
    assert og_results.active(cfg(base))
    assert og_results.run_root("s3://b/irv/", "r-r1") == "s3://b/irv/results/open-gira/r-r1"
    assert og_results.parse_filters("osm_extract=planet-latest; slice=1") == {
        "osm_extract": "planet-latest", "slice": "1"}
    assert og_results.parse_filters(float("nan")) == {}
    assert og_results.fetched_path(cfg(base), "transport_rail_edges", True) == (
        f"raw_data/open-gira/{fx.RUN}/transport_rail_edges.geoparquet")


def test_open_run_checks_the_run(base, tmp_path):
    cat = og_results.open_run(cfg(base))
    assert cat.version("rail_edges") == fx.RUN
    with pytest.raises(og_results.ResultsError, match="no published run"):
        og_results.open_run(cfg(base, run="20990101-0000000-r1"))


def test_fetch_layer_directory_and_single(base, tmp_path):
    cat = og_results.open_run(cfg(base))
    paths = og_results.fetch_layer(cat, "road_edges_damages", tmp_path / "roads")
    assert sorted(p.name for p in paths) == [
        f"{fx.EXTRACT}_aqueduct-river_slice-0.geoparquet",
        f"{fx.EXTRACT}_aqueduct-river_slice-1.geoparquet",
    ]
    (rail,) = og_results.fetch_layer(cat, "rail_edges", tmp_path / "rail.geoparquet",
                                     "osm_extract=planet-latest", single=True)
    gdf = gpd.read_parquet(rail)
    assert gdf.crs.name == "WGS 84" and list(gdf.edge_id) == ["rail_0", "rail_1"]
    with pytest.raises(og_results.ResultsError, match="3 items match"):
        og_results.fetch_layer(cat, "admin_boundaries", tmp_path / "x.geoparquet", single=True)
    with pytest.raises(og_results.ResultsError, match="no items match"):
        og_results.fetch_layer(cat, "rail_edges", tmp_path / "y", "osm_extract=europe-latest")


def test_results_lock(base, tmp_path):
    lock = og_results.lock(cfg(base), {"rail_edges": "", "admin_boundaries": "admin_level=2"},
                           tmp_path / "og.lock.yaml")
    from irv_datapkg.catalog import LockFile, check
    import yaml

    data = yaml.safe_load(lock.read_text())
    assert data["run"] == fx.RUN and "release" not in data
    assert list(data["datasets"]["admin_boundaries"]["items"]) == ["admin-level-2"]
    # valid against irv-datapkg's lock schema where it knows results locks
    from irv_datapkg.catalog import validate
    if "run_id" in validate.schema()["$defs"]:
        assert check(data, "lock") == []
    assert LockFile.load(lock).version("rail_edges") == fx.RUN

    locked = og_results.results_config({"open_gira_results": {"lock": str(lock)}})
    assert locked["run"] == fx.RUN and og_results.active(locked)
    cat = og_results.open_run(locked)
    (path,) = og_results.fetch_layer(cat, "admin_boundaries", tmp_path / "adm2.geoparquet", single=True)
    assert "GID_2" in gpd.read_parquet(path).columns
    with pytest.raises(og_results.ResultsError, match="pins run"):
        og_results.open_run({**locked, "run": "other-r1"})
    # a lock refuses changed content
    data["datasets"]["rail_edges"]["items"][fx.EXTRACT]["data"]["sha256"] = "0" * 64
    lock.write_text(yaml.safe_dump(data))
    from irv_datapkg.catalog import ChecksumError
    with pytest.raises(ChecksumError):
        og_results.fetch_layer(og_results.open_run(locked), "rail_edges", tmp_path / "z.gpq", single=True)


def test_cli_fetch(base, tmp_path):
    out = tmp_path / "nodes.geoparquet"
    assert og_results.main(["fetch", "--base", str(base), "--run", fx.RUN, "--layer", "rail_nodes",
                            "--single", "--output", str(out)]) == 0
    assert set(gpd.read_parquet(out).tag_railway) == {"station", "halt"}


# -- nbs-adaptation ------------------------------------------------------------


def notebook_columns():
    """Hazard columns listed in the opportunities notebook's Extract cell."""
    path = ETL / "pipelines/nbs-adaptation/extract-nbs-adaptation-opportunities.ipynb"
    source = next(
        "".join(c["source"]) for c in json.loads(path.read_text())["cells"]
        if "".join(c["source"]).startswith("read = [")
    )
    out = {}
    for block in re.split(r'\(\s*"', source)[1:]:
        name = block.split('"', 1)[0]
        out[name] = re.findall(r"'(hazard-[^']+)'", block)
    return out


def test_expected_columns_match_notebook():
    listed = notebook_columns()
    for kind in nbs.NBS_TYPES:
        assert sorted(nbs.expected_columns(kind)) == sorted(listed[f"{kind}_with_EAD"]), kind


def test_boundary_geojsonseq_and_nbs(base, tmp_path):
    cat = og_results.open_run(cfg(base))
    (hb,) = og_results.fetch_layer(cat, "hydrobasins_gadm", tmp_path / "hb.geoparquet", single=True)
    polys, points = tmp_path / "hb.geojsonld", tmp_path / "hb_points.geojsonld"
    nbs.boundary_geojsonseq(hb, "hybas_lev12_v1c_with_gadm_codes", polys, points)
    first = json.loads(polys.read_text().splitlines()[0])
    assert first["properties"]["HYBAS_ID"] == 11
    assert first["properties"]["bbox_wkt"].startswith("POLYGON")
    assert json.loads(points.read_text().splitlines()[0])["geometry"]["type"] == "Point"

    (opp,) = og_results.fetch_layer(cat, "nbs_opportunities", tmp_path / "m.gpq",
                                    {"nbs": "mangrove"}, single=True)
    out = nbs.nbs_with_ead([opp], "mangrove", tmp_path / "mangrove_with_EAD.geoparquet")
    gdf = gpd.read_parquet(out)
    assert set(nbs.expected_columns("mangrove")) <= set(gdf.columns)
    assert "bbox" not in gdf.columns


# -- Snakefiles ----------------------------------------------------------------


def legacy_results(root: Path) -> Path:
    """An open-gira results directory as the notebooks read it."""
    import geopandas
    from shapely.geometry import box

    for n in range(3):
        cols = {f"GID_{i}": ["A"] for i in range(n + 1)} | {f"NAME_{i}": ["a"] for i in range(n + 1)}
        path = root / f"input/admin-boundaries/admin-level-{n}.geoparquet"
        path.parent.mkdir(parents=True, exist_ok=True)
        geopandas.GeoDataFrame(cols, geometry=[box(0, 0, 1, 1)], crs="EPSG:4326").to_parquet(path)
    hb = geopandas.GeoDataFrame(
        {"HYBAS_ID": [7], **{f"{k}_{i}": ["A"] for i in range(3) for k in ("GID", "NAME")}},
        geometry=[box(0, 0, 1, 1)], crs="EPSG:4326").set_index("HYBAS_ID")
    path = root / "input/hydrobasins/hybas_lev12_v1c_with_gadm_codes.geoparquet"
    path.parent.mkdir(parents=True, exist_ok=True)
    hb.to_parquet(path)
    for kind in nbs.NBS_TYPES:
        for s, rows in ((0, 1), (1, 0)):  # slice-1 is empty, without HYBAS_ID
            path = root / f"slices/planet-latest_nbs/slice-{s}/{kind}_with_EAD.parquet"
            path.parent.mkdir(parents=True, exist_ok=True)
            cols = {"HYBAS_ID": [7] * rows, "area_ha": [1.0] * rows} if rows else {"x": []}
            geopandas.GeoDataFrame(cols, geometry=[box(0, 0, 1, 1)] * rows, crs="EPSG:4326").to_parquet(path)
    return root


@pytest.fixture
def workdir(tmp_path):
    """etl/-like working directory: pipelines/ and scripts/ linked from the repo."""
    for name in ("pipelines", "scripts"):
        os.symlink(ETL / name, tmp_path / name)
    return tmp_path


def run(pipeline, workdir, *args, check=True, **config):
    cfg_args = []
    for key, value in config.items():
        cfg_args.append(f"{key}={json.dumps(value) if isinstance(value, dict) else value}")
    return snakemake(
        ["-s", str(ETL / f"pipelines/{pipeline}/Snakefile"), "-d", str(workdir), "--cores", "1",
         *args, *(["--config", *cfg_args] if cfg_args else [])],
        check=check,
    )


def network_config(pipeline, csv="network_layers.csv"):
    return {"network_layers": str(ETL / f"pipelines/{pipeline}/{csv}"),
            "network_tilelayers": str(ETL / f"pipelines/{pipeline}/network_tilelayers.csv")}


@pytest.mark.parametrize("pipeline,csv,tiles", [
    ("osm_roads", "network_layers.csv", "network_tilelayers.csv"),
    ("osm_rail", "network_layers_edges.csv", "network_tilelayers_edges.csv"),
    ("osm_rail", "network_layers_nodes.csv", "network_tilelayers_nodes.csv"),
])
def test_network_pipelines_read_results_run(base, workdir, pipeline, csv, tiles):
    layers = {"network_layers": str(ETL / f"pipelines/{pipeline}/{csv}"),
              "network_tilelayers": str(ETL / f"pipelines/{pipeline}/{tiles}")}
    n_tiles = len(open(layers["network_tilelayers"]).read().strip().splitlines()) - 1
    # legacy: no fetch, networks_to_db reads the configured path
    out = run(pipeline, workdir, "-n", "-p", "--forceall", **layers).stdout
    assert "fetch_og_results" not in out
    assert jobs(out)["networks_to_db"] == n_tiles
    # results run: one fetch per network layer, feeding networks_to_db
    og = {"open_gira_results": {"base": str(base), "run": fx.RUN}}
    out = run(pipeline, workdir, "-n", "--forceall", **layers, **og).stdout
    assert jobs(out)["fetch_og_results"] == 1
    assert jobs(out)["networks_to_db"] == n_tiles
    # and the fetch itself runs
    import pandas
    ref = pandas.read_csv(layers["network_layers"]).ref[0]
    target = og_results.fetched_path(og_results.results_config(og), ref, pipeline == "osm_rail")
    run(pipeline, workdir, target, **layers, **og)
    fetched = workdir / target
    files = sorted(fetched.glob("*.geoparquet")) if fetched.is_dir() else [fetched]
    assert files and all(gpd.read_parquet(f).crs.name == "WGS 84" for f in files)
    if pipeline == "osm_roads":
        assert len(files) == 2


def test_nbs_adaptation_from_results_run(base, workdir):
    og = {"open_gira_results": {"base": str(base), "run": fx.RUN}}
    out = run("nbs-adaptation", workdir, "-n", **og).stdout
    assert jobs(out) == Counter({"boundary": 4, "boundary_geojsonseq": 4, "boundary_mbtiles": 8,
                                 "nbs_opportunities": 3, "all": 1})
    targets = ["raw_data/boundaries/admin-level-1_points.geojsonld",
               "raw_data/boundaries/hybas_lev12_v1c_with_gadm_codes.geojsonld",
               "raw_data/nbs-adaptation/mangrove_with_EAD.geoparquet"]
    run("nbs-adaptation", workdir, *targets, **og)
    gdf = gpd.read_parquet(workdir / targets[2])
    assert set(nbs.expected_columns("mangrove")) <= set(gdf.columns)


def test_nbs_adaptation_legacy_fallback(workdir, tmp_path):
    legacy = legacy_results(tmp_path / "og-results")
    targets = ["raw_data/boundaries/hybas_lev12_v1c_with_gadm_codes_points.geojsonld",
               "raw_data/nbs-adaptation/river_basin_afforestation_with_EAD.geoparquet"]
    run("nbs-adaptation", workdir, *targets, legacy_results_dir=str(legacy))
    first = json.loads((workdir / targets[0]).read_text().splitlines()[0])
    assert first["properties"]["HYBAS_ID"] == 7
    gdf = gpd.read_parquet(workdir / targets[1])
    assert len(gdf) == 1 and set(nbs.expected_columns("river_basin_afforestation")) <= set(gdf.columns)
