"""Load network features from a single source file-layer to database.
"""

import pandas
from sqlalchemy import delete
from sqlalchemy.orm import Session
from sqlalchemy.sql import exists
import geopandas

from backend.db.database import SessionLocal
from backend.db.models import Feature, FeatureLayer


def yield_features(layer, network_tile_layer):
    """Read from geoparquet into DB"""
    # Could limit columns here if we didnt need them all in JSON
    filters = [(layer.asset_type_column, "=", network_tile_layer.asset_type)]
    print(f"Loading geoparquet with filters: {filters}")
    gdf = geopandas.read_parquet(layer.path, filters=filters)
    print(f"Selected {gdf.shape} features from input")
    if gdf.crs.name != "WGS 84":
        raise Exception("Input geoparquet must be in WGS 84")
    gdf.reset_index()
    gdf["wkt"] = gdf.geometry.to_wkt()

    rename = {
        layer.asset_id_column: "asset_id",
        layer.asset_type_column: "asset_type",
    }
    remove = ["geometry"]

    print("Yielding Features...")
    for _, feature in gdf.iterrows():
        props = clean_props(feature, rename, remove=remove)
        props["sector"] = network_tile_layer.sector
        props["subsector"] = network_tile_layer.subsector

        yield Feature(
            string_id=props["asset_id"],
            layer=network_tile_layer.layer,
            properties=props,
            geom=feature.wkt,
        )


def clean_props(props: pandas.Series, rename: dict, remove=[]):
    clean = {}
    props = props.fillna("")
    props = props.to_dict()
    for k, v in props.items():
        if k in rename:
            clean[rename[k]] = v
        elif k in rename.values():
            clean[f"_{k}"] = v
        elif k in remove:
            continue

        else:
            clean[k] = v
    return clean


def get_network_layer(layer_name, network_layers):
    try:
        return network_layers[network_layers.ref == layer_name].iloc[0]
    except IndexError as e:
        print(f"Could not find {layer_name} in network layers.")
        raise e


def get_network_layer_by_ref(
    network_tile_layer_ref: str, network_layers: pandas.DataFrame
):
    try:
        return network_layers[network_layers.ref == network_tile_layer_ref].iloc[0]
    except IndexError as e:
        print(f"Could not find {network_tile_layer_ref} in network layers.")
        raise e


def get_network_layer_path(layer):
    return f"{layer.path}"


def get_tilelayer_by_layer_ref(layer_ref: str, network_tilelayers: pandas.DataFrame):
    return network_tilelayers[network_tilelayers.ref == layer_ref].iloc[0]


def get_tilelayer_by_layer_name(layer_name: str, network_tilelayers: pandas.DataFrame):
    return network_tilelayers[network_tilelayers.layer == layer_name].iloc[0]


def load_tile_feature_layer(db: Session, network_tile_layer):
    """Load FeatureLayer to DB if it doesnt exist"""
    layer_exists = db.query(
        exists().where(FeatureLayer.layer_name == network_tile_layer.layer)
    ).scalar()
    if not layer_exists:
        feature_layer = FeatureLayer(
            layer_name=network_tile_layer.layer,
            sector=network_tile_layer.sector,
            subsector=network_tile_layer.subsector,
            asset_type=network_tile_layer.asset_type,
        )
        db.add(feature_layer)
        db.commit()


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

    print("Layer", layer)
    network_tile_layer = get_tilelayer_by_layer_name(layer, network_tilelayers)
    print("Network TileLayer:", network_tile_layer)
    network_layer = get_network_layer_by_ref(network_tile_layer.ref, network_layers)
    print("Network Layer", network_layer)

    db: Session
    with SessionLocal() as db:
        print("Adding FeatureLayer if required")
        load_tile_feature_layer(db, network_tile_layer)
        print("Removing existing features for this tilelayer")
        db.execute(delete(Feature).where(Feature.layer == network_tile_layer.layer))
        db.commit()

        print("BEGINNING INGEST...")
        for i, feature in enumerate(yield_features(network_layer, network_tile_layer)):
            db.add(feature)
            if i % 1000 == 0:
                print("Dumped rows to DB, done:", i)
                db.commit()
        db.commit()
        print("Final Commit")

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{network_layer.path}\n\n")
        fh.write(f"Details:\n{str(network_layer)}\n")
