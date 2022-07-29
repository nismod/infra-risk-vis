"""Load network features from a single source file-layer to database.
"""
import fiona
import pandas

from pyproj import CRS, Transformer
from shapely.geometry import shape
from shapely.ops import transform
from sqlalchemy import delete
from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import Feature


def yield_features(layer, network_tilelayer, analysis_data_dir):
    """Read from file layer to modelled Features"""
    with fiona.open(get_network_layer_path(layer, analysis_data_dir), layer=layer.gpkg_layer) as src:
        from_crs = src.crs
        to_crs = CRS.from_epsg(4326)
        t = Transformer.from_crs(from_crs, to_crs, always_xy=True).transform

        def to_2d(x, y, z=None):
            return (x, y)

        rename = {
            layer.asset_id_column: "asset_id",
            layer.asset_type_column: "asset_type",
            layer.asset_min_cost_column: "cost_min",
            layer.asset_max_cost_column: "cost_max",
            layer.asset_mean_cost_column: "cost_mean",
            layer.asset_cost_unit_column: "cost_unit",
        }

        for i, feature in enumerate(src):
            if "id" in feature and "fid" not in feature["properties"]:
                feature["properties"]["fid"] = int(feature["id"])
            geom = transform(t, shape(feature["geometry"]))
            if geom.has_z:
                geom = transform(to_2d, geom)
            props = clean_props(feature["properties"], rename)
            # FIXME in the data
            if layer.ref == "transport_rail_edges":
                props["asset_type"] = "track"
            if layer.ref == "energy_edges":
                props["asset_type"] = "transmission"
            if layer.ref == "energy_targets":
                props["asset_type"] = "demand"
            props["asset_type"] = props["asset_type"].lower()
            props["asset_id"] = str(props["asset_id"])
            tilelayer_details = get_tilelayer_by_asset_type(layer.ref, props, network_tilelayers)
            props["sector"] = tilelayer_details.sector
            props["subsector"] = tilelayer_details.subsector

            yield Feature(
                id=int(layer.base_id + i),
                string_id=f"{props['asset_id']}",
                layer=tilelayer_details.layer,
                properties=props,
                geom=geom.wkt,
            )


def clean_props(props, rename):
    clean = {}
    for k, v in props.items():
        if k in rename:
            clean[rename[k]] = v
        elif k in rename.values():
            clean[f"_{k}"] = v
        else:
            clean[k] = v
    return clean


def get_network_layer(layer_name, network_layers):
    try:
        return network_layers[network_layers.ref == layer_name].iloc[0]
    except IndexError as e:
        print(f"Could not find {layer_name} in network layers.")
        raise e


def get_network_layer_path(layer, analysis_data_dir):
    return f"{analysis_data_dir}/{layer.path}"


def get_tilelayer_by_asset_type(layer_ref, props, network_tilelayers):
    try:
        return network_tilelayers[
            (network_tilelayers.asset_type == props["asset_type"])
            & (network_tilelayers.ref == layer_ref)
        ].iloc[0]
    except IndexError as e:
        print(f"Could not find {props['asset_type']} in tilelayers.")
        raise e


def get_tilelayers_by_network_source(network_source, network_tilelayers):
    try:
        return network_tilelayers[network_tilelayers.ref == network_source.ref].layer
    except IndexError as e:
        print(f"Could not find {network_source} in tilelayers.")
        raise e


if __name__ == "__main__":
    try:
        layer = snakemake.wildcards.layer
        output = snakemake.output
        analysis_data_dir = snakemake.config["analysis_data_dir"]

        network_layers = pandas.read_csv(snakemake.config["network_layers"])
        network_tilelayers = pandas.read_csv(snakemake.config["network_tilelayers"])

    except NameError:
        print("Expected to run from snakemake")
        exit()

    layer = get_network_layer(layer, network_layers)
    tilelayers = list(get_tilelayers_by_network_source(layer, network_tilelayers))

    db: Session
    with SessionLocal() as db:
        db.execute(delete(Feature).where(Feature.layer.in_(tilelayers)))
        db.commit()

        for i, feature in tqdm(enumerate(yield_features(layer, network_tilelayers, analysis_data_dir)), total=layer["count"]):
            db.add(feature)
            if i % 1000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(
            f"From:\n{get_network_layer_path(layer, analysis_data_dir)}|{layer.gpkg_layer}\n\n"
        )
        fh.write(f"Details:\n{str(layer)}\n")
