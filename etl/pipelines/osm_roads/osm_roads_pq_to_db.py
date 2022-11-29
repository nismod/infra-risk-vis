"""Load network features from a single source file-layer to database.
"""
import os
import sys
from typing import List
import warnings
from datetime import datetime

warnings.filterwarnings("error")

import pandas as pd
from sqlalchemy import delete, insert
from sqlalchemy.orm import Session
from sqlalchemy.sql import exists
import geopandas
from pyarrow import dataset
import pyarrow.dataset as ds
from shapely import wkb

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
    pq_fpath: str,
    filter: dataset.Expression,
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

    feature_columns = ["string_id", "layer", "properties", "geom"]

    rename = {
        layer.asset_id_column: "asset_id",
        layer.asset_type_column: "asset_type",
    }
    # Will remove columns from props containing these names
    remove = ["geometry", "hazard"]

    # Load PQ
    dataset = ds.dataset(pq_fpath, format="parquet")
    # Skip if null slice
    if dataset.count_rows() == 0:
        print(f"Skipping slice {pq_fpath} - zero rows (before filter)")
        return 0, 0, 0

    for idx, batch in enumerate(
        dataset.to_batches(batch_size=batch_size, filter=filter)
    ):
        #print(
        #    f"{datetime.now().isoformat()} - Doing batch {idx}, total features so far: {total_features}"
        #)
        if batch.num_rows == 0:
            continue
        try:
            # Load data to memory from disk
            df = batch.to_pandas()
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
                    },
                ),
                axis=1,
            )
            # Generate the string_id column
            df["string_id"] = df[layer.asset_id_column]
            df["layer"] = network_tile_layer.layer
            df["geom"] = df["geometry"].apply(lambda x: wkb_to_wkt(x))
            # Generate dicts for bulk insert
            feature_rows = df[feature_columns].to_dict(orient="records")
            db: Session
            with SessionLocal() as db:
                # This return_defaults makes it slower - but I cant see another way around it
                db.bulk_insert_mappings(Feature, feature_rows, return_defaults=True)
                db.commit()
            # Get the pks
            feature_pk_ids = [row["id"] for row in feature_rows]
            # Check we inserted enough rows
            if not len(feature_pk_ids) == batch.num_rows:
                print(f" Warning - Missed some rows in batch: {idx}")
            total_features += len(feature_pk_ids)

            # Generate a DF to join of the pk ids.  The list sorting should be enough, but we'll join for safety
            df_pk_ids = pd.DataFrame(
                [
                    {"id": item["id"], layer.asset_id_column: item["string_id"]}
                    for item in feature_rows
                ]
            )
            df_pk_ids.set_index(layer.asset_id_column, inplace=True)
            df_pks_joined = df.join(df_pk_ids)
            #print(
            #    "Joined PK Shapes:",
            #    df.shape[0],
            #    len(feature_rows),
            #    df_pks_joined.shape[0],
            #)

            # Generate RPDamages rows for this batch
            df_rp_damage = parse_rp_damage_batch(df_pks_joined, primary_key_column="id")
            if df_rp_damage.shape[0] > 0:
                # Impute a rename cols as required
                df_rp_damage.rename(
                    columns={
                        "MIN": "damage_amin",
                        "MEAN": "damage_mean",
                        "MAX": "damage_amax",
                    },
                    inplace=True,
                )
                df_rp_damage.reset_index(inplace=True)
                df_rp_damage.rename(columns={"id": "feature_id"}, inplace=True)
                df_rp_damage["loss_amin"] = 0.0
                df_rp_damage["loss_mean"] = 0.0
                df_rp_damage["loss_amax"] = 0.0
                df_rp_damage["exposure"] = 0.0
                # print("RP Damages Shape:", df_rp_damage.shape)
                rp_damage_rows = df_rp_damage.to_dict("records")
                db: Session
                with SessionLocal() as db:
                    # Dont need ids here
                    db.bulk_insert_mappings(
                        ReturnPeriodDamage, rp_damage_rows, return_defaults=False
                    )
                    db.commit()
                total_rp_damages += len(rp_damage_rows)
                del rp_damage_rows
                del df_rp_damage
            # Generate Expected Damages rows for this batch
            df_expected_damage = parse_exp_damage_batch(
                df_pks_joined, primary_key_column="id"
            )
            if df_expected_damage.shape[0] > 0:
                df_expected_damage.reset_index(inplace=True)
                df_expected_damage.rename(columns={"id": "feature_id"}, inplace=True)
                # print("Expected Damages Shape:", df_expected_damage.shape)
                expected_damage_rows = df_expected_damage.to_dict("records")
                db: Session
                with SessionLocal() as db:
                    # Dont need ids here
                    db.bulk_insert_mappings(
                        ExpectedDamage, expected_damage_rows, return_defaults=False
                    )
                    db.commit()
                    total_expected_damages += len(expected_damage_rows)
                del expected_damage_rows
                del df_expected_damage
        except Exception as err:
            print(f"Failed batch {idx} if pq file {pq_fpath}, due to {err} - skipped")

    return total_features, total_rp_damages, total_expected_damages


def parse_rp_damage_batch(df: pd.DataFrame, primary_key_column="id") -> pd.DataFrame:
    """Parse pandas df for rp damages"""

    data_cols = [c for c in df.columns if "rp" in c]

    melted = (
        df.melt(id_vars=primary_key_column, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    meta = melted.variable.str.extract(r"^hazard-(\w+)_(\w+)_(\w+)_(\d+)_(\w+)?")
    meta.columns = ["hazard", "rcp", "var", "epoch", "rp"]
    meta["var"].fillna("none", inplace=True)
    meta.loc[meta["hazard"] == "inunriver", "hazard"] = "river"
    meta.loc[meta["hazard"] == "inuncoast", "hazard"] = "coastal"
    meta["hazard"] = meta["hazard"].str.slice(0, 8)
    meta.loc[meta["rcp"] == "historical", "rcp"] = "baseline"
    meta["epoch"] = meta["epoch"].astype(int)
    meta["rp"] = meta["rp"].str.replace("rp", "").astype(int)

    return (
        melted.join(meta)
        .drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch", "rp"],
            columns="var",
            values="value",
        )
        .fillna(0)
    )


def parse_exp_damage_batch(df: pd.DataFrame, primary_key_column="id") -> pd.DataFrame:
    """Parse pandas df for rp damages"""
    data_cols = [c for c in df.columns if "EAD" in c]

    data = (
        df.melt(id_vars=primary_key_column, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    # parse string key column for metadata
    meta = data.variable.str.extract(r"^hazard-(\w+)_(\w+)_(\w+)_(\d+)")
    meta.columns = ["hazard", "rcp", "var", "epoch"]
    meta.loc[meta["hazard"] == "inunriver", "hazard"] = "river"
    meta.loc[meta["hazard"] == "inuncoast", "hazard"] = "coastal"
    meta["hazard"] = meta["hazard"].str.slice(0, 8)
    meta.loc[meta["rcp"] == "historical", "rcp"] = "baseline"

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
        columns={"MIN": "ead_amin", "MEAN": "ead_mean", "MAX": "ead_amax"}, inplace=True
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
    with open(network_layer.layer+'.txt', "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{network_layer.path}\n\n")
        fh.write(
            f"Details:\n{str(network_layer)}, total_feats: {total_feats}, total_rp_damages:{total_rp_damages}, total_expected_damages:{total_expected_damages}\n"
        )


class Layer:
    asset_id_column = "edge_id"
    asset_type_column = "tag_highway"
    layer = "road_edges_metro"


class NetworkTileLayer:
    asset_type = "secondary"
    sector = "transport"
    subsector = "road"
    layer = "road_edges_metro"


if __name__ == "__main__":
    # For testing
    layer = Layer()
    network_tile_layer = NetworkTileLayer()

    #try:
    #    layer = snakemake.wildcards.layer
    #    output = snakemake.output
    #    analysis_data_dir = snakemake.config["analysis_data_dir"]

    #    network_layers = pd.read_csv(snakemake.config["network_layers"])
    #    network_tilelayers = pd.read_csv(snakemake.config["network_tilelayers"])

    #th open("/opt/infra-risk/class_c_completed.txt", "r") as f:except NameError:
    #    print("Expected to run from snakemake")
    #    exit()

    #print("Layer", layer)
    #network_tile_layer = get_tilelayer_by_layer_name(layer, network_tilelayers)
    #print("Network TileLayer:", network_tile_layer)
    #network_layer = get_network_layer_by_ref(network_tile_layer.ref, network_layers)
    #print("Network Layer", network_layer)
    #print(f"PQ Path: {network_layer.path}")

    #filter = ds.field(network_layer.asset_type_column) == network_tile_layer.asset_type
    filter = ds.field(layer.asset_type_column) == network_tile_layer.asset_type

    #pq_fpath = "/home/dusted/code/oxford/infra-risk-vis/etl/raw_data/processed_data/input/osm_roads/20221102_egypt_with_EAD_and_RP_test.geoparquet"
    pq_dir = "/opt/infra-risk/etl/raw_data/processed_data/input/osm_roads/20221103_global_road_EAD_and_cost_per_RP"

    db: Session
    with SessionLocal() as db:
        print("Adding FeatureLayer if required")
        load_tile_feature_layer(db, network_tile_layer)

    with open("/opt/infra-risk/class_metro_completed.txt", "r") as f:
        skip = set([line.strip() for line in f])

    dryrun = False

    # Walk dir and load slices
    total_features = 0
    total_rp_damages = 0
    total_expected_damages = 0
    for root, dirs, files in os.walk(pq_dir, topdown=False):
        for file in files:
            if '.geoparquet' not in file:
                print ('skipping non parquet file: ', file)
                continue

            if file in skip:
                print('skipping already completed: ', file)
                continue
            else:
                if dryrun:
                    print('would not skip: ', file)

            if dryrun:
                continue

            pq_fpath = os.path.join(root, file)
            print(
                f"{datetime.now().isoformat()} - Processing {os.path.basename(pq_fpath)}"
            )

            file_features, file_rp_damages, file_expected_damages = load_batches(
                layer,
                network_tile_layer,
                pq_fpath,
                filter,
                batch_size=10,
            )
            total_features += file_features
            total_rp_damages += file_rp_damages
            total_expected_damages += file_expected_damages
            with open("/opt/infra-risk/class_metro_completed.txt", "w") as f:
                f.write(f'{file}\n')
            print(
                f"{datetime.now().isoformat()} - Completed {os.path.basename(pq_fpath)}, total features: {total_features}, total_rp_damages: {total_rp_damages}, total_expected_damages: {total_expected_damages}"
            )

    print(f"Done - Loaded {total_features} features")
    write_output_log(
        layer, total_features, total_rp_damages, total_expected_damages
    )
