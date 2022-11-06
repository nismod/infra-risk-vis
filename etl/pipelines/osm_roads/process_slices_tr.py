import json
import logging
import os
from glob import glob

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from tqdm import tqdm

pq_dir = '/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/20221103_global_road_EAD_and_cost_per_RP/'

# # Count rows
# fnames = glob(f"{pq_dir}/*.geoparquet")
# for f in fnames:
#     pf = pq.ParquetFile(f)
#     print(os.path.basename(f), pf.metadata.num_rows)

# # Fails with schema mismatch
# dataset = ds.dataset(pq_dir, format="parquet")
# dataset = pq.ParquetDataset(pq_dir)
# breakpoint()



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

    data = (
        melted.join(meta)
        .drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch", "rp"],
            columns="var",
            values="value",
        )
        .fillna(0)
        .reset_index()
        .rename(
            columns={
                "id": "feature_id",
                "MIN": "damage_amin",
                "MEAN": "damage_mean",
                "MAX": "damage_amax",
            },
        )
    )
    # ensure all columns are present - may be missing in case the data didn't
    # have any non-zero values in this batch
    expected_columns = [
        "damage_amin",
        "damage_mean",
        "damage_amax",
        "loss_amin",
        "loss_mean",
        "loss_amax",
        "exposure",
    ]
    ensure_columns(data, expected_columns)
    ordered_cols = [
        "feature_id", "hazard", "rcp", "epoch", "rp", "exposure",
        "damage_amin", "damage_mean", "damage_amax",
        "loss_amin", "loss_mean", "loss_amax",
    ]
    return data[ordered_cols]



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

    # pivot back up so we end with a row per uid, hazard etc. (see index columns below)
    # and columns for each damage type, each with min/mean/max
    data = (
        data.drop(columns="variable")
        .pivot(
            index=[primary_key_column, "hazard", "rcp", "epoch"],
            columns=["var"],
            values="value",
        )
        .fillna(0)
        .reset_index()
        .rename(
            columns={
                "id": "feature_id",
                "MIN": "ead_amin",
                "MEAN": "ead_mean",
                "MAX": "ead_amax",
            },
        )
    )

    # ensure all columns are present - may be missing in case the data didn't
    # have any non-zero values in this batch
    expected_columns = [
        "protection_standard",
        "ead_amin",
        "ead_mean",
        "ead_amax",
        "eael_amin",
        "eael_mean",
        "eael_amax",
    ]
    ensure_columns(data, expected_columns)
    ordered_cols = [
        "feature_id", "hazard", "rcp", "epoch", "protection_standard",
        "ead_amin", "ead_mean", "ead_amax",
        "eael_amin", "eael_mean", "eael_amax",
    ]
    return data[ordered_cols]



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
    return json.dumps(res)

# base_id = 0
# # total is 1_708_288
# class Layer:
#     asset_id_column = "edge_id"
#     asset_type_column = "tag_highway"
#     layer = "road_edges_motorway"

# class NetworkTileLayer:
#     asset_type = "motorway"
#     sector = "transport"
#     subsector = "road"
#     layer = "road_edges_motorway"

base_id = 2_000_000
# total is 2_598_000
class Layer:
    asset_id_column = "edge_id"
    asset_type_column = "tag_highway"
    layer = "road_edges_trunk"

class NetworkTileLayer:
    asset_type = "trunk"
    sector = "transport"
    subsector = "road"
    layer = "road_edges_trunk"

# base_id = 5_000_000
# # total is 4_727_438
# class Layer:
#     asset_id_column = "edge_id"
#     asset_type_column = "tag_highway"
#     layer = "road_edges_primary"

# class NetworkTileLayer:
#     asset_type = "primary"
#     sector = "transport"
#     subsector = "road"
#     layer = "road_edges_primary"

# base_id =  10_000_000
# # total is  6_447_961
# class Layer:
#     asset_id_column = "edge_id"
#     asset_type_column = "tag_highway"
#     layer = "road_edges_secondary"

# class NetworkTileLayer:
#     asset_type = "secondary"
#     sector = "transport"
#     subsector = "road"
#     layer = "road_edges_secondary"

# base_id = 20_000_000
# # total is 9_251_742
# class Layer:
#     asset_id_column = "edge_id"
#     asset_type_column = "tag_highway"
#     layer = "road_edges_tertiary"

# class NetworkTileLayer:
#     asset_type = "tertiary"
#     sector = "transport"
#     subsector = "road"
#     layer = "road_edges_tertiary"


if __name__ == '__main__':
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    layer = Layer()
    network_tile_layer = NetworkTileLayer()

    filter = ds.field(layer.asset_type_column) == network_tile_layer.asset_type
    batch_size=100_000

    fnames = list(sorted(glob(f"{pq_dir}/*.geoparquet")))
    # test with largest:
    # fnames = list(sorted(glob(f"{pq_dir}/*538.geoparquet")))

    for f in tqdm(fnames):
        dataset = ds.dataset(f, format="parquet")

        feature_columns = ["id", "string_id", "layer", "properties", "geom"]

        rename = {
            layer.asset_id_column: "asset_id",
            layer.asset_type_column: "asset_type",
        }
        # Will remove columns from props containing these names
        remove = ["geometry", "hazard"]

        for idx, batch in enumerate(
            dataset.to_batches(batch_size=batch_size, filter=filter)
        ):
            if batch.num_rows == 0:
                continue

            df = batch.to_pandas()
            df["id"] = np.arange(base_id, base_id + batch.num_rows)
            base_id += batch.num_rows
            logging.info(f"Max id: {base_id}, current batch: {batch.num_rows}")

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
            # Write out
            df["string_id"] = df[layer.asset_id_column]
            df["layer"] = "road_edges_" + df[layer.asset_type_column]
            csv_path = f.replace(".geoparquet", f"_{network_tile_layer.layer}_features_{idx}.csv")
            df["geom"] = ""
            df[feature_columns].to_csv(csv_path, index=False)

            df_rp_damage = parse_rp_damage_batch(df, primary_key_column="id")
            rp_csv_path = f.replace(".geoparquet", f"_{network_tile_layer.layer}_rp_{idx}.csv")
            df_rp_damage.to_csv(rp_csv_path, index=False)

            df_expected_damage = parse_exp_damage_batch(df, primary_key_column="id")
            exp_csv_path = f.replace(".geoparquet", f"_{network_tile_layer.layer}_exp_{idx}.csv")
            df_expected_damage.to_csv(exp_csv_path, index=False)

        logging.info(f"Completed {f}:{network_tile_layer.layer}")
