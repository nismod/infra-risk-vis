import ujson as json
import logging

import numpy as np
import pyarrow.dataset as ds
from pygeos.io import from_wkb, to_geojson


def feature_as_geojson(row, hazard_cols):
    properties = {
        "asset_id": row["asset_id"],
        "asset_type": row["asset_type"],
    }

    for hazard_col in hazard_cols:
        # e.g.     'hazard-inunriver_rcp4p5_MAX_2080_EAD'
        # pattern: "^hazard-(\w+)_(\w+)_(\w+)_(\d+)_EAD"
        #              ["hazard", "rcp", "var", "epoch"]
        hazard, rcp, _, epoch, _ = hazard_col.split("_")
        hazard = hazard.replace("hazard-inun", "") # hardcoded "inun"!
        epoch = int(epoch)
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



fname = "/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/20221104_global_rail_EAD_and_cost_per_RP.geoparquet"
filter=None
base_id = 60_000_000
class Layer:
    asset_id_column = "edge_id"
    asset_type_column = "tag_railway"
    layer = "rail_edges"

class NetworkTileLayer:
    asset_type = "rail"
    sector = "transport"
    subsector = "rail"
    layer = "rail_edges"


# fname = "/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/20221027_global_rail_stations.geoparquet"
# base_id = 65_000_000
# filter = ds.field("tag_railway") == "station"
# class Layer:
#     asset_id_column = "node_id"
#     asset_type_column = "tag_railway"
#     layer = "rail_nodes"

# class NetworkTileLayer:
#     asset_type = "station"
#     sector = "transport"
#     subsector = "rail"
#     layer = "rail_nodes"


if __name__ == '__main__':
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)

    batch_size=100_000

    layer = Layer()
    network_tile_layer = NetworkTileLayer()

    asset_id_column = layer.asset_id_column
    output_dir = '/mnt/d/mert2014/Users/mert2014/data/infra-tmp/'

    with open(f"{output_dir}/{network_tile_layer.layer}.geojson", "w") as fh:
        dataset = ds.dataset(fname, format="parquet")
        feature_columns = ["id", "asset_id", "asset_type", "geometry"]
        rename = {
            asset_id_column: "asset_id",
        }

        for idx, batch in enumerate(
            dataset.to_batches(batch_size=batch_size, filter=filter)
        ):
            if batch.num_rows == 0:
                continue
            df = batch.to_pandas().reset_index().rename(rename, axis='columns')
            df["asset_type"] = network_tile_layer.layer

            # pattern: "^hazard-(\w+)_(\d+)_(\w+)_(\w+)_ead"
            #              ["hazard", "epoch", "rcp", "var"]
            hazard_columns = [c for c in df.columns if "EAD" in c and "MAX" in c]

            df["id"] = np.arange(base_id, base_id + batch.num_rows)
            df["geometry"] = to_geojson(from_wkb(df.geometry))
            base_id += batch.num_rows
            logging.info(f"Max id: {base_id}, current batch: {batch.num_rows}")
            process_df(df[feature_columns+hazard_columns], fh, hazard_columns)

        logging.info(f"Completed {fname}:{network_tile_layer.layer}")
