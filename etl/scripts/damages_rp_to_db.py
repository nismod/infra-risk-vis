"""Load return period damages from a single file to database
"""
import pyarrow.parquet as pq

from sqlalchemy.orm import Session
from tqdm import tqdm

from backend.db.database import SessionLocal
from backend.db.models import ReturnPeriodDamage


def yield_return_period_damages(fname):
    """Generate ReturnPeriodDamage database model objects, reading through
    potentially very large parquet files in row batches.
    """
    batch_size = 100
    rp_pf = pq.ParquetFile(fname)
    rp_batches = rp_pf.iter_batches(batch_size)

    for batch in rp_batches:
        batch_df = batch.to_pandas()
        # in case of data not having non-zero values in this batch
        expected_columns = [
            "damage_amin",
            "damage_mean",
            "damage_amax",
            # "loss_amin",
            # "loss_mean",
            # "loss_amax",
            "exposure",
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
                loss_amin=0,
                loss_mean=0,
                loss_amax=0,
            )


def ensure_columns(data, expected_columns):
    for col in expected_columns:
        if col not in data.columns:
            data[col] = 0
    return data


if __name__ == "__main__":
    try:
        input = snakemake.input
        fname = snakemake.input.fname
        layer = snakemake.wildcards.layer
        output = snakemake.output
    except NameError:
        print("Expected to run from snakemake")
        exit()

    damage_rps = yield_return_period_damages(fname)

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
