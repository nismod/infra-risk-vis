import ujson as json
import logging
import os
from glob import glob

import numpy as np
import pandas as pd
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from pygeos.io import from_wkb, to_geojson
# from shapely.geometry import mapping
from tqdm import tqdm

pq_dir = '/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/20221103_global_road_EAD_and_cost_per_RP/'


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
        "ead_amin",
        "ead_mean",
        "ead_amax",
    ]
    ensure_columns(data, expected_columns)
    ordered_cols = [
        "feature_id", "hazard", "rcp", "epoch",
        "ead_amin", "ead_mean", "ead_amax",
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

# base_id = 2_000_000
# # total is 2_598_000
# class Layer:
#     asset_id_column = "edge_id"
#     asset_type_column = "tag_highway"
#     layer = "road_edges_trunk"

# class NetworkTileLayer:
#     asset_type = "trunk"
#     sector = "transport"
#     subsector = "road"
#     layer = "road_edges_trunk"

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

base_id = 20_000_000
# total is 9_251_742
class Layer:
    asset_id_column = "edge_id"
    asset_type_column = "tag_highway"
    layer = "road_edges_tertiary"

class NetworkTileLayer:
    asset_type = "tertiary"
    sector = "transport"
    subsector = "road"
    layer = "road_edges_tertiary"


def feature_as_geojson(row, hazard_cols):
    properties = {
        "asset_id": row["asset_id"],
        "asset_type": row["asset_type"],
    }

    for hazard_col in hazard_cols:
        hazard, rcp, _, epoch, _ = hazard_col.split("_")
        hazard = hazard.replace("hazard-inun", "") # hardcoded "inun"!
        key = f"{hazard}__rcp_{rcp}__epoch_{epoch}"
        properties[f"ead__{key}"] = row[hazard_col]

    return {
        "type": "Feature",
        "id": row["id"],
        "geometry": json.loads(row["geometry"]),
        "properties": properties,
    }


def process_df(df, fh, hazard_columns):
    for row in df.to_dict(orient="records"):
        feature = feature_as_geojson(row, hazard_columns)
        json.dump(feature, fh, indent=0, ensure_ascii=False)
        fh.write("\n")


if __name__ == '__main__':
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    layer = Layer()
    network_tile_layer = NetworkTileLayer()

    filter = ds.field(layer.asset_type_column) == network_tile_layer.asset_type
    batch_size=100_000

    fnames = list(sorted(glob(f"{pq_dir}/*.geoparquet")))
    # test with largest:
    # fnames = list(sorted(glob(f"{pq_dir}/*538.geoparquet")))

    output_dir = '/mnt/d/mert2014/Users/mert2014/data/infra-tmp/'
    with open(f"{output_dir}/{network_tile_layer.layer}.geojson", "w") as fh:
        for f in fnames:
            dataset = ds.dataset(f, format="parquet")
            feature_columns = ["id", "asset_id", "asset_type", "geometry"]
            rename = {
                "asset_type": "paved",
                layer.asset_id_column: "asset_id",
                layer.asset_type_column: "asset_type",
            }

            for idx, batch in enumerate(
                dataset.to_batches(batch_size=batch_size, filter=filter)
            ):
                if batch.num_rows == 0:
                    continue
                df = batch.to_pandas().rename(rename, axis='columns')

                # pattern: "^hazard-(\w+)_(\w+)_(\w+)_(\d+)_EAD"
                #              ["hazard", "rcp", "var", "epoch"]
                hazard_columns = [c for c in df.columns if "EAD" in c and "MAX" in c]

                df["id"] = np.arange(base_id, base_id + batch.num_rows)
                df["geometry"] = to_geojson(from_wkb(df.geometry))
                base_id += batch.num_rows
                logging.info(f"Max id: {base_id}, current batch: {batch.num_rows}")
                process_df(df[feature_columns+hazard_columns], fh, hazard_columns)

            logging.info(f"Completed {f}:{network_tile_layer.layer}")
