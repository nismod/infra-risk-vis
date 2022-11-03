"""Load network features from a single source file-layer to database.
"""
import os
import sys
from typing import List

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
) -> None:
    """
    Load the given parquet file into the Db in batches

        # Read parquet schema using dataset

        # Iterate batches of rows to memory

        # Generate dictionaries for each feature row in the batch
        # Load the features and get the ids (converting wkb with shapely if req) using insert of dict mappings

        # Generate dictionaries for damages and expected rows in batch (multiple rows per DF row)

        # Load Damages and Expected
    """
    total_batches = 0
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

    for idx, batch in enumerate(
        dataset.to_batches(batch_size=batch_size, filter=filter)
    ):
        print(f"Doing batch {idx}")
        if batch.num_rows == 0:
            continue
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
        df_pks_joined = df.join(df_pk_ids, rsuffix="_pk")
        print(
            "Joined PK Shapes:", df.shape[0], len(feature_rows), df_pks_joined.shape[0]
        )
        df_pks_joined.to_pickle(
            "/home/dusted/code/oxford/infra-risk-vis/etl/pipelines/osm_roads/df_pks_joined.pkl"
        )

        # Generate RPDamages rows for this batch
        df_rp_damage = parse_rp_damage_batch(df_pks_joined, primary_key_column="id")
        print("RP Damages Shape:", df_rp_damage.shape)
        expected_damage_rows = df_rp_damage.to_dict(orient="records")
        print(df_rp_damage)
        if df_rp_damage.shape[0] > 0:
            break
        # db: Session
        # with SessionLocal() as db:
        #     # Dont need ids here
        #     db.bulk_insert_mappings(
        #         ReturnPeriodDamage, rp_damage_rows, return_defaults=False
        #     )
        #     db.commit()

        # # Generate Expected Damages rows for this batch
        # df_expected_damage = parse_exp_damage_batch(
        #     df_pks_joined, primary_key_column="id"
        # )
        # print("Expected Damages Shape:", df_expected_damage.shape)
        # expected_damage_rows = df_expected_damage.to_dict(orient="records")
        # db: Session
        # with SessionLocal() as db:
        #     # Dont need ids here
        #     db.bulk_insert_mappings(
        #         ExpectedDamage, expected_damage_rows, return_defaults=False
        #     )
        #     db.commit()

    print(total_features)


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

    return (
        melted.join(meta)
        .drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch", "rp"],
            columns="var",
            values="value",
        )
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

    data.rename({"min": "ead_amin", "mean": "ead_mean", "max": "ead_amax"})
    data["eael_amin"] = 0
    data["eael_mean"] = 0
    data["eael_max"] = 0

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
                remove_substr is True
                and any([item.find(k) != -1 for item in remove]) is False
            ):
                continue
            else:
                if k in remove:
                    continue
            clean[k] = v
    # Merge additional keys
    # res = {**clean, **append}
    # print(res)
    return clean


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


class Layer:
    asset_id_column = "edge_id"
    asset_type_column = "tag_highway"
    layer = "road_edges_motorway"


class NetworkTileLayer:
    asset_type = "motorway"
    sector = "transport"
    subsector = "road"
    layer = "road_edges_motorway"


if __name__ == "__main__":
    # layer = pd.DataFrame(
    #     [
    #         {
    #             "asset_id_column": "edge_id",
    #             "asset_type_column": "tag_highway",
    #             "layer": "road_edges_motorway",
    #         }
    #     ]
    # )

    # network_tile_layer = pd.DataFrame(
    #     [
    #         {
    #             "asset_type": "motorway",
    #             "sector": "transport",
    #             "subsector": "road",
    #             "layer": "road_edges_motorway",
    #         }
    #     ]
    # )

    layer = Layer()
    network_tile_layer = NetworkTileLayer()

    filter = ds.field(layer.asset_type_column) == network_tile_layer.asset_type

    pq_fpath = "/home/dusted/code/oxford/infra-risk-vis/etl/raw_data/processed_data/input/osm_roads/20221102_egypt_with_EAD_and_RP_test.geoparquet"

    db: Session
    with SessionLocal() as db:
        print("Adding FeatureLayer if required")
        load_tile_feature_layer(db, network_tile_layer)
        print("Removing existing features for this tilelayer")
        db.execute(delete(Feature).where(Feature.layer == network_tile_layer.layer))
        db.commit()

    load_batches(
        layer,
        network_tile_layer,
        pq_fpath,
        filter,
        batch_size=1000,
    )
