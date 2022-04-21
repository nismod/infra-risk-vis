"""Load return period damages from a single file to database
"""
import pyarrow.parquet as pq

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import ReturnPeriodDamage


def yield_return_period_damages(exposure_fname, damage_fname, loss_fname):
    """Generate ReturnPeriodDamage database model objects, reading through
    potentially very large parquet files in row batches.
    """
    batch_size = 100
    exp_pf = pq.ParquetFile(exposure_fname)
    exp_batches = exp_pf.iter_batches(batch_size)
    dmg_pf = pq.ParquetFile(damage_fname)
    dmg_batches = dmg_pf.iter_batches(batch_size)
    loss_pf = pq.ParquetFile(loss_fname)
    loss_batches = loss_pf.iter_batches(batch_size)

    # Loop while there are still batches, joining each set of three batches
    # (one from each file) at each iteration. This works so long as the Parquet
    # files share asset order and number of rows.
    while True:
        try:
            exp_df = parse_rp_damage_batch(next(exp_batches)).rename(
                columns={"none": "exposure"}
            )
            dmg_df = parse_rp_damage_batch(next(dmg_batches)).rename(
                columns={
                    "amin": "damage_amin",
                    "mean": "damage_mean",
                    "amax": "damage_amax",
                }
            )
            loss_df = parse_rp_damage_batch(next(loss_batches)).rename(
                columns={
                    "amin": "loss_amin",
                    "mean": "loss_mean",
                    "amax": "loss_amax",
                }
            )
            batch_df = exp_df.join(dmg_df).join(loss_df).fillna(0).reset_index()

            # in case of data not having non-zero values in this batch
            expected_columns = [
                'damage_amin', 'damage_mean', 'damage_amax',
                'loss_amin', 'loss_mean', 'loss_amax',
                'exposure'
            ]
            ensure_columns(batch_df, expected_columns)

            for row in batch_df.itertuples():
                yield ReturnPeriodDamage(
                    feature_id=row.uid,
                    hazard=row.hazard,
                    rcp=row.rcp,
                    epoch=row.epoch,
                    rp=row.rp,
                    exposure=row.exposure,
                    damage_amin=row.damage_amin,
                    damage_mean=row.damage_mean,
                    damage_amax=row.damage_amax,
                    loss_amin=row.loss_amin,
                    loss_mean=row.loss_mean,
                    loss_amax=row.loss_amax,
                )
        except StopIteration:
            break


def parse_rp_damage_batch(batch):
    """Parse a parquet (arrow) row batch to pandas
    """
    data = batch.to_pandas()
    data_cols = [c for c in batch.schema.names if "rp" in c]

    melted = (
        data.melt(id_vars="uid", value_vars=data_cols)
        .query("value > 0")
        .reset_index(drop=True)
    )

    meta = melted.variable.str.extract(
        r"^(\w+)__rp_(\d+)__rcp_([\w\d.]+)__epoch_(\d+)__?conf_([^_]+)_?(\w+)?"
    )
    meta.columns = ["hazard", "rp", "rcp", "epoch", "conf", "var"]
    meta["var"].fillna("none", inplace=True)

    return (
        melted.join(meta)
        .drop(columns="variable")
        .pivot(
            index=["uid", "hazard", "rp", "rcp", "epoch", "conf"],
            columns="var",
            values="value",
        )
        .query("conf == '50' or conf == 'None'")
    )


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

    damage_rps = yield_return_period_damages(
        exposure, damage, loss
    )

    db: Session
    with SessionLocal() as db:
        for i, damage_rp in enumerate(
            tqdm(damage_rps, desc=f"{layer}_rp")
        ):
            db.add(damage_rp)
            if i % 10000 == 0:
                db.commit()
        db.commit()

    with open(str(output), "w") as fh:
        fh.write(f"Loaded to database.\n\n")
        fh.write(f"From:\n{input}\n\n")
