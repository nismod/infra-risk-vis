"""Load network features from a single source file-layer to database.
"""
import logging
import os
import sys
import warnings
from datetime import datetime

warnings.filterwarnings("error")

import numpy as np
import pandas as pd
from sqlalchemy.orm import Session
from sqlalchemy.sql import exists
from pyarrow import dataset
import pyarrow.dataset as ds
from shapely import wkb
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from common.db.database import SessionLocal
from common.db.models import Feature, FeatureLayer, ReturnPeriodDamage, ExpectedDamage


def load_pq(pq_fpath: str) -> dataset:
    """
    Load PQ File to data set (doesnt load to memory)
    """
    return ds.dataset(pq_fpath, format="parquet")


def pq_equality_expression(column: str, value: str) -> dataset.Expression:
    """
    Generate a pyarrow equality expression for the given column and value
    """
    return ds.field(column) == value


def wkb_to_wkt(geom_wkb) -> str:
    return wkb.loads(geom_wkb).wkt


def load_batches(
    layer,
    network_tile_layer,
    base_id,
    pq_fpath: str,
    batch_size: int = 100,
) -> int:
    """
    Load the given parquet file into the Db in batches

        # Read parquet schema using dataset

        # Iterate batches of rows to memory

        # Generate dictionaries for each feature row in the batch
        # Load the features and get the ids (converting wkb with shapely if req) using insert of dict mappings

        # Generate dictionaries for damages and expected rows in batch (multiple rows per DF row)

        # Load Damages and Expected
    """
    total_features = 0
    total_rp_damages = 0
    total_expected_damages = 0

    feature_columns = ["id", "string_id", "layer", "properties"] # , "geom"

    rename = {
        layer.asset_id_column: "asset_id",
    }
    # Will remove columns from props containing these names
    remove = ["geometry", "hazard"]

    # Load PQ
    dataset = ds.dataset(pq_fpath, format="parquet")
    # Skip if null slice
    if dataset.count_rows() == 0:
        print(f"Skipping slice {pq_fpath} - zero rows")
        return 0, 0, 0

    for idx, batch in tqdm(enumerate(
        dataset.to_batches(batch_size=batch_size)
    ), total=(3678243//batch_size)):
        # print(
        #     f"{datetime.now().isoformat()} - Doing batch {idx}, total features so far: {total_features}"
        # )
        if batch.num_rows == 0:
            continue

        # Load data to memory from disk
        df = batch.to_pandas().reset_index()
        df["id"] = np.arange(base_id, base_id + batch.num_rows)
        base_id += batch.num_rows
        # logging.info(f"Max id: {base_id}, current batch: {batch.num_rows}")
        # Map the feature rows to as-expected by the table
        # - Generate the props field
        df["properties"] = df.apply(
            lambda x: clean_props(
                x,
                rename,
                remove,
                remove_substr=True,
                append={
                    "sector": network_tile_layer.sector,
                    "subsector": network_tile_layer.subsector,
                    "asset_type": network_tile_layer.asset_type,
                },
            ),
            axis=1,
        )
        # Generate the string_id column
        df["string_id"] = df[layer.asset_id_column]
        df["layer"] = network_tile_layer.layer
        # df["geom"] = df["geometry"].apply(lambda x: wkb_to_wkt(x))

        # Generate dicts for bulk insert
        feature_rows = df[feature_columns].to_dict(orient="records")
        total_features += len(feature_rows)

        db: Session
        with SessionLocal() as db:
            db.bulk_insert_mappings(Feature, feature_rows, return_defaults=False)
            db.commit()

        # Generate RPDamages rows for this batch
        df_rp_damage = parse_rp_damage_batch(df, primary_key_column="id")
        if df_rp_damage.shape[0] > 0:
            df_rp_damage.reset_index(inplace=True)
            df_rp_damage.rename(columns={"id": "feature_id"}, inplace=True)
            df_rp_damage["loss_amin"] = 0.0
            df_rp_damage["loss_mean"] = 0.0
            df_rp_damage["loss_amax"] = 0.0
            df_rp_damage["exposure"] = 0.0

            rp_damage_rows = df_rp_damage.to_dict(orient="records")
            db: Session
            with SessionLocal() as db:
                # Dont need ids here
                db.bulk_insert_mappings(
                    ReturnPeriodDamage, rp_damage_rows, return_defaults=False
                )
                db.commit()
            total_rp_damages += len(rp_damage_rows)

        # Generate Expected Damages rows for this batch
        df_expected_damage = parse_exp_damage_batch(
            df, primary_key_column="id"
        )
        if df_expected_damage.shape[0] > 0:
            df_expected_damage.reset_index(inplace=True)
            df_expected_damage.rename(columns={"id": "feature_id"}, inplace=True)
            # print("Expected Damages Shape:", df_expected_damage.shape)
            expected_damage_rows = df_expected_damage.to_dict(orient="records")
            db: Session
            with SessionLocal() as db:
                # Dont need ids here
                db.bulk_insert_mappings(
                    ExpectedDamage, expected_damage_rows, return_defaults=False
                )
                db.commit()
                total_expected_damages += len(expected_damage_rows)

    return total_features, total_rp_damages, total_expected_damages


def parse_rp_damage_batch(df: pd.DataFrame, primary_key_column="id") -> pd.DataFrame:
    """Parse pandas df for rp damages"""

    data_cols = [c for c in df.columns if "rp" in c]

    melted = (
        df.melt(id_vars=primary_key_column, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    meta = melted.variable.str.extract(r"^hazard-(\w+)_(\d+)_(\w+)_rp(\d+)_(\w+)")
    meta.columns = ["hazard", "epoch", "rcp", "rp", "var"]
    meta["var"].fillna("none", inplace=True)
    meta["hazard"] = "cyclone"
    meta["epoch"] = meta["epoch"].astype(int)
    meta.loc[meta["rcp"] == "rcp8p5", "rcp"] = "8.5"
    meta["rp"] = meta["rp"].astype(int)

    data = (
        melted.join(meta)
        .drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch", "rp"],
            columns="var",
            values="value",
        )
        .fillna(0)
    )

    data.rename(
        columns={"min": "damage_amin", "mean": "damage_mean", "max": "damage_amax"}, inplace=True
    )
    return data


def parse_exp_damage_batch(df: pd.DataFrame, primary_key_column="id") -> pd.DataFrame:
    """Parse pandas df for rp damages"""
    data_cols = [c for c in df.columns if "ead" in c]

    data = (
        df.melt(id_vars=primary_key_column, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    # parse string key column for metadata
    meta = data.variable.str.extract(r"^hazard-(\w+)_(\d+)_(\w+)_(\w+)_ead")
    meta.columns = ["hazard", "epoch", "rcp", "var"]
    meta["hazard"] = "cyclone"
    meta["epoch"] = meta["epoch"].astype(int)
    meta.loc[meta["rcp"] == "rcp8p5", "rcp"] = "8.5"

    # join metadata columns
    data = data.join(meta)
    if "protection_standard" not in df.columns:
        data["protection_standard"] = 0
    else:
        data.loc[data.defended == "undefended", "protection_standard"] = 0

    # pivot back up so we end with a row per uid, hazard etc. (see index columns below)
    # and columns for each damage type, each with min/mean/max
    data = (
        data.drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch", "protection_standard"],
            columns=["var"],
            values="value",
        )
        .fillna(0)
    )

    data.rename(
        columns={"min": "ead_amin", "mean": "ead_mean", "max": "ead_amax"}, inplace=True
    )

    # ensure all columns are present - may be missing in case the data didn't
    # have any non-zero values in this batch
    expected_columns = [
        "ead_amin",
        "ead_mean",
        "ead_amax",
        "eael_amin",
        "eael_mean",
        "eael_amax",
    ]
    ensure_columns(data, expected_columns)

    return data.reset_index()


def ensure_columns(data, expected_columns):
    for col in expected_columns:
        if col not in data.columns:
            data[col] = 0
    return data


def clean_props(props, rename: dict, remove=[], remove_substr=True, append: dict = {}):
    """
    Generate cleaned (renamed / removed properties)
    ::kwarg substr bool If True then any column with the strings
        from remove list as substring will be removed
    ::kwarg append dict Optionall merge some additional properties at the end
    """
    clean = {}
    props = props.fillna("")
    props = props.to_dict()
    for k, v in props.items():
        if k in rename:
            clean[rename[k]] = v
        elif k in rename.values():
            clean[f"_{k}"] = v
        else:
            if (
                (remove_substr is True)
                and (any([k.find(item) != -1 for item in remove]))
            ):
                continue
            else:
                if k in remove:
                    continue
            clean[k] = v
    # Merge additional keys
    res = {**clean, **append}
    return res


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


def get_network_layer(layer_name, network_layers):
    try:
        return network_layers[network_layers.ref == layer_name].iloc[0]
    except IndexError as e:
        print(f"Could not find {layer_name} in network layers.")
        raise e


def get_network_layer_by_ref(network_tile_layer_ref: str, network_layers: pd.DataFrame):
    try:
        return network_layers[network_layers.ref == network_tile_layer_ref].iloc[0]
    except IndexError as e:
        print(f"Could not find {network_tile_layer_ref} in network layers.")
        raise e


def get_network_layer_path(layer):
    return f"{layer.path}"


def get_tilelayer_by_layer_ref(layer_ref: str, network_tilelayers: pd.DataFrame):
    return network_tilelayers[network_tilelayers.ref == layer_ref].iloc[0]


def get_tilelayer_by_layer_name(layer_name: str, network_tilelayers: pd.DataFrame):
    return network_tilelayers[network_tilelayers.layer == layer_name].iloc[0]


def write_output_log(
    network_layer, total_feats: int, total_rp_damages: int, total_expected_damages: int
):
    with open(str(f"{network_layer}"), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{network_layer.path}\n\n")
        fh.write(
            f"Details:\n{str(network_layer)}, total_feats: {total_feats}, total_rp_damages:{total_rp_damages}, total_expected_damages:{total_expected_damages}\n"
        )


base_id = 50_000_000
class Layer:
    asset_id_column = "edge_id"
    asset_type_column = ""
    layer = "power_transmission"


class NetworkTileLayer:
    asset_type = "transmission"
    sector = "power"
    subsector = "transmission"
    layer = "power_transmission"


if __name__ == "__main__":
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    # For testing
    network_layer = Layer()
    network_tile_layer = NetworkTileLayer()

    pq_fpath = "/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/grid_damages.geoparquet"
    network_layer.path = pq_fpath

    db: Session
    with SessionLocal() as db:
        print("Adding FeatureLayer if required")
        load_tile_feature_layer(db, network_tile_layer)

    # Walk dir and load slices
    total_features = 0
    total_rp_damages = 0
    total_expected_damages = 0

    print(
        f"{datetime.now().isoformat()} - Processing {os.path.basename(pq_fpath)}"
    )

    file_features, file_rp_damages, file_expected_damages = load_batches(
        network_layer,
        network_tile_layer,
        base_id,
        pq_fpath,
        batch_size=10000,
    )
    total_features += file_features
    total_rp_damages += file_rp_damages
    total_expected_damages += file_expected_damages
    print(
        f"{datetime.now().isoformat()} - Completed {os.path.basename(pq_fpath)}, total features: {total_features}, total_rp_damages: {total_rp_damages}, total_expected_damages: {total_expected_damages}"
    )

    print(f"Done - Loaded {total_features} features")
    write_output_log(
        network_layer, total_features, total_rp_damages, total_expected_damages
    )
