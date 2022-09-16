"""Load return period damages from a single file to database
"""
import pandas
import pyarrow.parquet as pq

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import ReturnPeriodDamage


def yield_return_period_damages(exposure_fname, damage_fname, loss_fname):
    """Generate ReturnPeriodDamage database model objects, reading through
    potentially very large parquet files in row batches.
    """
    exp_df = pandas.read_parquet(exposure_fname)
    dmg_df = pandas.read_parquet(damage_fname)
    loss_df = pandas.read_parquet(loss_fname)

    exp_df = parse_rp_df(exp_df).rename(
        columns={"exposure_amax": "exposure"}
    )
    dmg_df = parse_rp_df(dmg_df).rename(
        columns={
            "min": "damage_amin",
            "mean": "damage_mean",
            "max": "damage_amax",
        }
    )
    loss_df = parse_rp_df(loss_df).rename(
        columns={
            "min": "loss_amin",
            "mean": "loss_mean",
            "max": "loss_amax",
        }
    )
    batch_df = exp_df.join(dmg_df).join(loss_df).fillna(0).reset_index()

    # in case of data not having non-zero values in this batch
    expected_columns = [
        "damage_amin",
        "damage_mean",
        "damage_amax",
        "loss_amin",
        "loss_mean",
        "loss_amax",
        "exposure",
    ]
    ensure_columns(batch_df, expected_columns)

    for row in batch_df.itertuples():
        if row.epoch == 1980:
            epoch = 2010
        else:
            epoch = row.epoch
        yield ReturnPeriodDamage(
            feature_id=row.uid,
            hazard=row.hazard,
            rcp=row.rcp,
            epoch=epoch,
            rp=row.rp,
            exposure=row.exposure,
            damage_amin=row.damage_amin,
            damage_mean=row.damage_mean,
            damage_amax=row.damage_amax,
            loss_amin=row.loss_amin,
            loss_mean=row.loss_mean,
            loss_amax=row.loss_amax,
        )


def parse_rp_df(data):
    """Parse a pandas dataframe to long format"""
    data_cols = [c for c in data.columns if "rp" in c]

    melted = (
        data.melt(id_vars="uid", value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    meta = melted.variable.str.extract(
        r"^(\w+)__rp_(\d+)__rcp_([\w\d.]+)__epoch_(\d+)__conf_([^_]+)__subs_([^_]+)__model_([^_]+)_?(\w+)?"
    )
    meta.columns = ["hazard", "rp", "rcp", "epoch", "conf", "subs", "model", "var"]

    long = (
        melted.join(meta)
        .drop(columns=["variable", "conf", "subs", "model", "var"])
        .groupby(["uid", "hazard", "rp", "rcp", "epoch"])
        .agg(['min', 'mean', 'max'])
        .fillna(0)
    )
    long.columns = [agg for (_, agg) in long.columns]
    return long



def ensure_columns(data, expected_columns):
    for col in expected_columns:
        if col not in data.columns:
            data[col] = 0
    return data


if __name__ == "__main__":
    try:
        input = snakemake.input
        exposure = snakemake.input.exposure
        damage = snakemake.input.damage
        loss = snakemake.input.loss
        layer = snakemake.wildcards.layer
        output = snakemake.output
    except NameError:
        print("Expected to run from snakemake")
        exit()

    damage_rps = yield_return_period_damages(exposure, damage, loss)

    db: Session
    with SessionLocal() as db:
        for i, damage_rp in enumerate(tqdm(damage_rps, desc=f"{layer}_rp")):
            db.add(damage_rp)
            if i % 10000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{input}\n\n")
