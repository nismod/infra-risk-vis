"""fetch_datapkg script against the synthetic catalogue, with and without a lock."""

import hashlib

import pytest
import yaml

import datapkg_fetch
import datapkg_mapping as dm
from irv_datapkg.catalog import Catalog, ChecksumError

KEY = "inunriver_rcp4p5_00000NorESM1-M_2030_rp00100"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_fetch_from_base_and_release(catalogue, etl_root, cache, tmp_path):
    out = tmp_path / "raster/raw/aqueduct" / f"{KEY}.tif"
    out.parent.mkdir(parents=True)
    item, path = datapkg_fetch.fetch_layer(
        "aqueduct", KEY, out, base=str(catalogue), release="2026.1", cache_dir=cache, root=etl_root)
    assert item.id == "inunriver_rcp4p5_noresm1-m_2030_rp00100"
    assert path == out and sha256(out) == item.asset.sha256
    assert not out.with_name(out.name + ".part").exists()


def test_fetch_with_lock(catalogue, etl_root, cache, tmp_path):
    lock = tmp_path / "datapkg.lock.yaml"
    Catalog.open(str(catalogue), release="2026.1").write_lock(lock, ["aqueduct_flood"])
    out = tmp_path / f"{KEY}.tif"
    main_args = ["--pipeline", "aqueduct", "--key", KEY, "--output", str(out),
                 "--lock", str(lock), "--cache-dir", cache]
    import os
    cwd = os.getcwd()
    os.chdir(etl_root)
    try:
        datapkg_fetch.main(main_args)
    finally:
        os.chdir(cwd)
    assert out.exists()

    # the lock's checksum is what the Snakemake param carries
    data = yaml.safe_load(lock.read_text())
    spec = {"dataset": "aqueduct_flood", "item": "inunriver_rcp4p5_noresm1-m_2030_rp00100"}
    assert dm.locked_sha256(data, spec) == sha256(out)

    # a tampered lock is refused by the client
    entry = data["datasets"]["aqueduct_flood"]["items"][spec["item"]]["data"]
    entry["sha256"] = "0" * 64
    lock.write_text(yaml.safe_dump(data))
    with pytest.raises(ChecksumError):
        datapkg_fetch.fetch_layer("aqueduct", KEY, out, lock=str(lock), cache_dir=cache,
                                  root=etl_root)


def test_spec_sha_mismatch_refused(catalogue, etl_root, cache, tmp_path):
    spec = {"dataset": "aqueduct_flood", "item": "inunriver_rcp4p5_noresm1-m_2030_rp00100",
            "sha256": "f" * 64}
    with pytest.raises(dm.MappingError, match="differs from the lock"):
        datapkg_fetch.fetch_layer("aqueduct", KEY, tmp_path / "x.tif", base=str(catalogue),
                                  cache_dir=cache, spec=spec, root=etl_root)
