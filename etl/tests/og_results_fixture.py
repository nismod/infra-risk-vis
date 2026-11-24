"""A tiny synthetic open-gira results run, {BASE}/results/open-gira/{run}/.

Written with the irv-datapkg producer API following irv-datapkg
docs/catalogue-spec.md section 8 (layer ids, dimensions, run root fields),
as open-gira's ``publish_results`` writes it. Nothing is downloaded.
"""

import json
from pathlib import Path

import geopandas as gpd
from shapely.geometry import LineString, Point, box

from irv_datapkg.catalog import CatalogWriter, build_collection, build_item

RUN = "20261124-0123456-r1"
COMMIT = "0123456789abcdef0123456789abcdef01234567"
EXTRACT = "planet-latest"
EAD = "hazard-inunriver_rcp4p5_MEAN_2030_EAD"
RP = "hazard-inunriver_rcp4p5_MEAN_2030_rp100"


def _gpq(gdf, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    gdf.to_parquet(path, index=False, schema_version="1.1.0", write_covering_bbox=True)
    return path


def _edges(prefix, tag, values):
    return gpd.GeoDataFrame(
        {"edge_id": [f"{prefix}_{i}" for i in range(len(values))], tag: values,
         EAD: [float(i) for i in range(len(values))], RP: [5.0] * len(values)},
        geometry=[LineString([(i, 0), (i + 0.5, 0.5)]) for i in range(len(values))],
        crs="EPSG:4326",
    )


def _boxes(n, **cols):
    return gpd.GeoDataFrame(cols, geometry=[box(i, 0, i + 1, 1) for i in range(n)], crs="EPSG:4326")


def layers(src: Path) -> dict:
    """{layer: (dimensions, [(item id, properties, file)])}."""
    admin = {}
    for n in range(3):
        cols = {}
        for i in range(n + 1):
            cols[f"GID_{i}"] = [f"AAA.{i}", f"BBB.{i}"]
            cols[f"NAME_{i}"] = [f"A{i}", f"B{i}"]
        admin[n] = _gpq(_boxes(2, **cols), src / f"admin-level-{n}.parquet")
    adm2 = {k: v for n in range(3) for k, v in {f"GID_{n}": ["AAA"], f"NAME_{n}": ["A"]}.items()}
    nbs = {}
    for name in ("river_basin_afforestation", "mangrove", "landslide_slope_vegetation"):
        nbs[name] = _gpq(
            _boxes(2, HYBAS_ID=[11, 12], GID_0=["AAA", "AAA"], GID_1=["a", "a"], GID_2=["b", "b"],
                   area_ha=[1.0, 2.0], **{f"hazard-x_{name}_EAD__road": [1.0, 0.0]}),
            src / f"nbs_{name}.parquet",
        )
    return {
        "road_edges_damages": (["osm_extract", "hazard", "slice"], [
            (f"{EXTRACT}_aqueduct-river_slice-{s}", {"osm_extract": EXTRACT, "hazard": "aqueduct-river", "slice": s},
             _gpq(_edges(f"road{s}", "tag_highway", ["motorway", "trunk", "primary"]), src / f"road{s}.parquet"))
            for s in (0, 1)
        ]),
        "rail_edges": (["osm_extract"], [
            (EXTRACT, {"osm_extract": EXTRACT}, _gpq(_edges("rail", "tag_railway", ["rail", "rail"]), src / "rail.parquet")),
        ]),
        "rail_nodes": (["osm_extract"], [
            (EXTRACT, {"osm_extract": EXTRACT}, _gpq(gpd.GeoDataFrame(
                {"node_id": ["n0", "n1"], "tag_railway": ["station", "halt"]},
                geometry=[Point(0, 0), Point(1, 1)], crs="EPSG:4326"), src / "nodes.parquet")),
        ]),
        "nbs_opportunities": (["osm_extract", "nbs"], [
            (f"{EXTRACT}_{name}", {"osm_extract": EXTRACT, "nbs": name}, path) for name, path in nbs.items()
        ]),
        "admin_boundaries": (["admin_level"], [
            (f"admin-level-{n}", {"admin_level": n}, path) for n, path in admin.items()
        ]),
        "hydrobasins_gadm": (["level"], [
            ("hybas_lev12_v1c_with_gadm_codes", {"level": 12},
             _gpq(_boxes(2, HYBAS_ID=[11, 12], **{k: v * 2 for k, v in adm2.items()}), src / "hybas.parquet")),
        ]),
    }


def build(root: Path, run: str = RUN) -> Path:
    """Write the run under ``{root}/base/results/open-gira/{run}``; return ``{root}/base``."""
    base = Path(root) / "base"
    run_root = base / "results" / "open-gira" / run
    writer = CatalogWriter(run_root)
    for layer, (dims, items) in layers(Path(root) / "src").items():
        writer.add_collection(build_collection(
            layer, title=layer, description=f"synthetic {layer}", license="ODbL-1.0",
            citation="synthetic", upstream_url="https://github.com/nismod/open-gira",
            upstream_version=COMMIT, dimensions=dims,
        ))
        for item_id, props, path in items:
            writer.add_item(build_item(
                path, layer, run, item_id,
                properties={**props, "irv:source_path": f"x/{path.name}"},
                source_url=f"open-gira://results/x/{path.name}", retrieved="2026-11-24", commit=COMMIT,
            ), path)
    doc = json.loads((run_root / "catalog.json").read_text())
    doc.update({
        "id": f"open-gira-{run}", "irv:producer": "open-gira", "irv:run": run,
        "irv:commit": COMMIT, "irv:dirty": False, "irv:config_hash": "ab" * 32,
        "irv:input_lock": None, "irv:created": "2026-11-24T12:00:00Z",
    })
    (run_root / "catalog.json").write_text(json.dumps(doc, indent=2))
    return base
