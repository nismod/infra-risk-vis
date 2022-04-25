"""Load expected damages from a single file to database
"""
import logging

import pandas

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import AdaptationCostBenefit


def yield_adaptation(data):
    """Generate NPVDamage database model objects."""
    for row in data.itertuples():
        yield AdaptationCostBenefit(
            feature_id=row.uid,
            hazard=row.hazard,
            rcp=row.rcp,
            adaptation_name=row.adaptation_option,
            adaptation_protection_level=row.protection_level,
            adaptation_cost=row.adapt_cost_npv,
            avoided_ead_amin=row.avoided_ead_amin,
            avoided_ead_mean=row.avoided_ead_mean,
            avoided_ead_amax=row.avoided_ead_amax,
            avoided_eael_amin=row.avoided_eael_amin,
            avoided_eael_mean=row.avoided_eael_mean,
            avoided_eael_amax=row.avoided_eael_amax,
        )


def parse_adaptation(data):
    """Parse to tidy (long) dataframe"""
    # files come in "fat" format, with string column names representing
    # hazard, RCP, variable (avoided_EAD/EAEL), and min/mean/max suffix
    data_cols = [c for c in data.columns if "EA" in c]

    data = data.rename(
        columns={
            "flood_depth_protection_level": "protection_level",
            "cyclone_damage_curve_reduction": "protection_level",
        }
    )
    # corner case for handling protection against "all" floods - set depth to 999
    if "flood_protection_level" in data.columns:
        data.loc[
            data.flood_protection_level == "All", "protection_level"
        ] = 999

    id_vars = ["uid", "adaptation_option", "protection_level", "adapt_cost_npv"]

    # melt to long format
    data = (
        data.melt(id_vars=id_vars, value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    # parse string key column for metadata
    meta = data.variable.str.extract(r"^([^_]+)__rcp_(\d.\d)__(\w+)$")
    meta.columns = ["hazard", "rcp", "var"]

    # join metadata columns
    data = data.join(meta)

    # pivot back up so we end with a row per uid, hazard etc. (see index columns below)
    # and columns for each damage type, each with min/mean/max
    data = (
        data.drop(columns="variable")
        .pivot(
            index=[
                "uid",
                "hazard",
                "rcp",
                "adaptation_option",
                "protection_level",
                "adapt_cost_npv",
            ],
            columns=["var"],
            values="value",
        )
        .fillna(0)
    )

    data.columns = [var.lower() for var in data.columns]

    # ensure all columns are present - may be missing in case the data didn't
    # have any non-zero values in this batch
    expected_columns = [
        "avoided_ead_amin",
        "avoided_ead_mean",
        "avoided_ead_amax",
        "avoided_eael_amin",
        "avoided_eael_mean",
        "avoided_eael_amax",
    ]
    ensure_columns(data, expected_columns)

    return data.reset_index()


def ensure_columns(data, expected_columns):
    for col in expected_columns:
        if col not in data.columns:
            logging.warning(f"Filling expected column '{col}' with zero")
            data[col] = 0
    return data


if __name__ == "__main__":
    try:
        input = snakemake.input
        avoided_risk_fname = snakemake.input.avoided_risk
        uid_fname = snakemake.input.uid
        layer_name = snakemake.wildcards.layer
        output = snakemake.output
        network_layers_fname = snakemake.config["network_layers"]
    except NameError:
        print("Expected to run from snakemake")
        exit()

    network_layers = pandas.read_csv(network_layers_fname)
    try:
        network_layer = network_layers[network_layers.damage_ref == layer_name].iloc[0]
    except IndexError as e:
        print(f"Could not find {layer_name} in network layers.")
        raise e

    adaptation = pandas.read_csv(avoided_risk_fname).set_index(
        network_layer.asset_id_column
    )
    ids = pandas.read_parquet(uid_fname).set_index(network_layer.asset_id_column)
    adaptation = adaptation.join(ids).reset_index(drop=True)

    adaptation_df = parse_adaptation(adaptation)
    adaptation = yield_adaptation(adaptation_df)

    db: Session
    with SessionLocal() as db:
        for i, damage_npv in enumerate(
            tqdm(adaptation, desc=f"{layer_name}_adaptation", total=len(adaptation_df))
        ):
            db.add(damage_npv)
            if i % 1000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{input}\n\n")
