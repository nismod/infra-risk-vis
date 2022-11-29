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
        # pattern: "^hazard-(\w+)_(\d+)_(\w+)_(\w+)_ead"
        #              ["hazard", "epoch", "rcp", "var"]
        hazard, epoch, rcp, var, _ = hazard_col.split("_")
        if var == "max":
            hazard = "cyclone" # hardcoded!
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


if __name__ == '__main__':
    logging.basicConfig(format="%(asctime)s %(message)s", level=logging.INFO)
    pq_fpath = "/home/tom/projects/infra-risk-vis/etl/raw_data/processed_data/input/grid_damages.geoparquet"

    network_tile_layers = [
        ("power_transmission", "openstreetmap"),
        ("power_distribution", "gridfinder"),
    ]
    asset_id_column = "edge_id"

    base_id = 50_000_000
    batch_size=100_000

    output_dir = '/mnt/d/mert2014/Users/mert2014/data/infra-tmp/'
    for network_tile_layer, source in network_tile_layers:
        filter = ds.field("source") == source

        with open(f"{output_dir}/{network_tile_layer}.geojson", "w") as fh:
            dataset = ds.dataset(pq_fpath, format="parquet")
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
                df["asset_type"] = network_tile_layer

                # pattern: "^hazard-(\w+)_(\d+)_(\w+)_(\w+)_ead"
                #              ["hazard", "epoch", "rcp", "var"]
                hazard_columns = [c for c in df.columns if "ead" in c and "max" in c]

                df["id"] = np.arange(base_id, base_id + batch.num_rows)
                df["geometry"] = to_geojson(from_wkb(df.geometry))
                base_id += batch.num_rows
                logging.info(f"Max id: {base_id}, current batch: {batch.num_rows}")
                process_df(df[feature_columns+hazard_columns], fh, hazard_columns)

            logging.info(f"Completed {pq_fpath}:{network_tile_layer}")
