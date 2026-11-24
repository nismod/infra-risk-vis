"""Workflow-level checks: rule precedence (dry runs) and a real fetch_datapkg run.

The dry runs use the real etl/ directory and touch nothing but .snakemake/.
The execution test runs the real Snakefile in a scratch working directory with
a 3-row aqueduct layers.csv and the synthetic catalogue: no database,
terracotta or network.
"""

import json
import os
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from conftest import ETL

pytest.importorskip("snakemake")

AQ = "raster/raw/aqueduct/inunriver_rcp4p5_00000NorESM1-M_2030_rp00100.tif"
STORM = "raster/raw/storm/STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_10_YR_RP.tif"
GHSL = "raster/raw/ghsl_pop/GHS_POP_E2020_GLOBE_R2023A_54009_1000_V1_0.tif"
JRC = "raster/raw/jrc_flood/floodMapGL_rp10y.tif"
WIRED = ["aqueduct", "storm", "ghsl_pop", "jrc_flood", "landslide_arup"]


def snakemake(args, cwd=ETL, check=True):
    env = dict(os.environ, PATH=f"{Path(sys.executable).parent}{os.pathsep}{os.environ['PATH']}")
    proc = subprocess.run(
        [sys.executable, "-m", "snakemake", "--nolock", *args],
        cwd=cwd, env=env, capture_output=True, text=True,
    )
    if check and proc.returncode != 0:
        raise AssertionError(proc.stdout + proc.stderr)
    return proc


def jobs(output):
    return Counter(re.findall(r"^rule (\S+):", output, flags=re.M))


def datapkg_config(**cfg):
    return ["--config", "datapkg=" + json.dumps(cfg)]


def test_disabled_by_default_uses_pipeline_rules():
    out = snakemake(["-n", AQ, STORM, GHSL]).stdout
    assert "fetch_datapkg" not in out
    assert jobs(out) == Counter({"aqueduct_download_raw_data": 1, "storm_download": 1,
                                 "storm_unpack": 1, "ghsl_pop_download_and_unzip_raw_data": 1})


def test_enabled_pipelines_fetch_from_catalogue():
    out = snakemake(["-n", AQ, STORM, GHSL, JRC,
                     *datapkg_config(enabled=WIRED, base="/nonexistent", release="2026.1")]).stdout
    assert jobs(out) == Counter({"fetch_datapkg": 4})


def test_switch_is_per_dataset():
    out = snakemake(["-n", AQ, STORM, *datapkg_config(enabled=["storm"], base="/x")]).stdout
    assert jobs(out) == Counter({"aqueduct_download_raw_data": 1, "fetch_datapkg": 1})


def test_full_dataset_dag_downstream_unchanged():
    out = snakemake(["-n", "raster/metadata/aqueduct.flag",
                     *datapkg_config(enabled=["aqueduct"], base="/x")]).stdout
    counts = jobs(out)
    assert counts["fetch_datapkg"] == counts["clip_raster"] == counts["cloud_optimise_raster"] == 349
    assert counts["ingest_rasters"] == counts["ingest_metadata"] == 1
    assert "aqueduct_download_raw_data" not in counts


def test_unlisted_raw_paths_are_not_claimed():
    # a raw path that is not a layers.csv row still goes to the pipeline rule
    target = "raster/raw/aqueduct/not_a_layer.tif"
    out = snakemake(["-n", target, *datapkg_config(enabled=["aqueduct"], base="/x")]).stdout
    assert jobs(out) == Counter({"aqueduct_download_raw_data": 1})


def test_fetch_datapkg_runs_against_local_catalogue(catalogue, etl_root, tmp_path):
    # --config must come last: it consumes the following arguments
    # scratch workdir: pipelines/ (subset layers), scripts/, ../envs/test/.etl.env
    os.symlink(ETL / "scripts", etl_root / "scripts")
    shutil.copy(ETL / "config.yml", etl_root / "config.yml")
    env_dir = etl_root.parent / "envs" / "test"
    env_dir.mkdir(parents=True)
    (env_dir / ".etl.env").write_text("IRV_DATAPKG_CACHE=" + str(tmp_path / "cache") + "\n")
    cfg = dict(enabled=["aqueduct"], base=str(catalogue), release="2026.1", lock="datapkg.lock.yaml")
    common = ["-s", str(ETL / "Snakefile"), "-d", str(etl_root), "--cores", "1",
              "--config", "environment=test", "datapkg=" + json.dumps(cfg)]

    # lock first (as committed in a real deployment), then fetch through it
    snakemake(["datapkg_lock", *common])
    lock = etl_root / "datapkg.lock.yaml"
    assert "aqueduct_flood" in lock.read_text()
    # record item ids in layers.csv from the locked catalogue
    import layers_from_catalogue
    assert layers_from_catalogue.main(["aqueduct", "--root", str(etl_root)]) == 0
    snakemake([AQ, *common])
    assert (etl_root / AQ).stat().st_size > 0

    # nothing to do on a second run
    out = snakemake(["-n", AQ, *common]).stdout
    assert jobs(out) == Counter()

    # with item ids recorded, the job carries the locked checksum (a re-pin
    # re-runs only the changed layers) and reads through the lock
    (etl_root / AQ).unlink()
    out = snakemake(["-n", "-p", AQ, *common]).stdout
    assert jobs(out) == Counter({"fetch_datapkg": 1})
    assert "'--lock' 'datapkg.lock.yaml'" in out and '"sha256": "' in out

    # CI check: the lock covers every layer with current checksums
    snakemake(["datapkg_lock_check", *common])
    # re-lock: with every row naming an item, only those items are locked
    lock.unlink()
    proc = snakemake(["-p", "datapkg_lock", *common])
    assert "'--item' 'aqueduct_flood/" in proc.stdout + proc.stderr
    import yaml
    locked = yaml.safe_load(lock.read_text())["datasets"]["aqueduct_flood"]["items"]
    assert len(locked) == 3
    snakemake(["datapkg_lock_check", *common])
    # a lock missing a needed item fails the check
    data = yaml.safe_load(lock.read_text())
    data["datasets"]["aqueduct_flood"]["items"].popitem()
    lock.write_text(yaml.safe_dump(data))
    assert snakemake(["datapkg_lock_check", *common], check=False).returncode != 0
