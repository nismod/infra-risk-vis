"""Source raw rasters from the irv-datapkg catalogue.

For every pipeline listed in ``datapkg.enabled`` (config.yml), ``fetch_datapkg``
writes ``raster/raw/{pipeline}/{key}.tif`` for each row of
``pipelines/{pipeline}/layers.csv``, from the catalogue item that row maps to
(``pipelines/{pipeline}/datapkg.yml``). Everything from ``clip_raster`` onwards
is unchanged.

The rule only claims the exact raw paths listed in the enabled pipelines'
layers.csv, and takes precedence (ruleorder) over any pipeline rule that writes
to ``raster/raw/{pipeline}/``. Pipelines not enabled keep their own download
rules. With nothing enabled, this file defines no rules.

N.B. include after the pipeline modules so their rules can be ordered.
"""

import json
import re
import sys

sys.path.insert(0, f"{workflow.basedir}/scripts")
import datapkg_mapping as datapkg

DATAPKG = datapkg.datapkg_config(config)
DATAPKG_LOCK = datapkg.read_lock(DATAPKG["lock"])


def _datapkg_specs():
    """{"pipeline/key": item spec} for every layer of every enabled pipeline."""
    specs = {}
    for pipeline in DATAPKG["enabled"]:
        mapping = datapkg.load_mapping(pipeline)
        layers = datapkg.read_layers(pipeline)
        for row in layers.to_dict(orient="records"):
            key = row["filename"][: -len(".tif")]
            spec = datapkg.item_spec(row, mapping)
            sha256 = datapkg.locked_sha256(DATAPKG_LOCK, spec)
            if sha256:
                spec["sha256"] = sha256  # a re-pinned checksum re-triggers the job
            specs[f"{pipeline}/{key}"] = spec
    return specs


DATAPKG_SPECS = _datapkg_specs() if DATAPKG["enabled"] else {}


def _datapkg_cli_args():
    """Catalogue options for scripts/datapkg_fetch.py.

    Both --base and --lock are passed when set: once the lock file exists,
    a configured ``base`` deliberately overrides the base recorded in the
    lock (irv-datapkg >= 0.3.0 ``Catalog.from_lock(base=...)``), so a
    deployment can read the same locked versions and checksums from a
    mirror. Leave ``base`` as recorded in the lock (or null) to use the
    lock's base; ``release`` is ignored while the lock exists.
    """
    args = []
    for option, key in (
        ("--base", "base"),
        ("--private-base", "private_base"),
        ("--release", "release"),
        ("--cache-dir", "cache_dir"),
    ):
        if DATAPKG.get(key):
            args += [option, str(DATAPKG[key])]
    if DATAPKG_LOCK is not None:
        args += ["--lock", DATAPKG["lock"]]
    return " ".join(f"'{a}'" for a in args)


def _datapkg_lock_args():
    """Selection for `irv-datapkg lock`: --item DATASET/ID for each layer that
    names an item; datasets with filter-resolved layers are locked whole."""
    whole, items = datapkg.lock_selection(DATAPKG_SPECS.values())
    args = [a for i in items for a in ("--item", i)] + whole
    return " ".join(f"'{a}'" for a in args)


if DATAPKG_SPECS:

    rule fetch_datapkg:
        """
        Fetch a raw raster from the irv-datapkg catalogue (sha256-verified).
        """
        output:
            "raster/raw/{DATAPKG_LAYER}.tif"
        wildcard_constraints:
            DATAPKG_LAYER="|".join(re.escape(k) for k in sorted(DATAPKG_SPECS)),
        params:
            spec=lambda wildcards: json.dumps(DATAPKG_SPECS[wildcards.DATAPKG_LAYER], sort_keys=True),
            pipeline=lambda wildcards: wildcards.DATAPKG_LAYER.split("/")[0],
            key=lambda wildcards: wildcards.DATAPKG_LAYER.split("/", 1)[1],
            catalogue=_datapkg_cli_args(),
        resources:
            disk_mb=1000,
        priority:
            70,
        shell:
            """
            python scripts/datapkg_fetch.py \
                --pipeline {params.pipeline} \
                --key {params.key} \
                --output {output} \
                --spec {params.spec:q} \
                {params.catalogue}
            """

    # take precedence over pipeline rules writing to raster/raw/{enabled}/
    _raw = re.compile(r"^raster/raw/([^/{}]+)/")
    for _rule in list(workflow.rules):
        if _rule.name == "fetch_datapkg":
            continue
        for _out in _rule.output:
            _m = _raw.match(str(_out))
            if _m and _m.group(1) in DATAPKG["enabled"]:
                workflow.ruleorder("fetch_datapkg", _rule.name)
                break


_DATAPKG_LOCK_PARAMS = dict(
    selection=_datapkg_lock_args(),
    base=DATAPKG["base"] or "",
    release=DATAPKG["release"] or "",
    private_base=f"--private-base '{DATAPKG['private_base']}'" if DATAPKG["private_base"] else "",
)


rule datapkg_lock:
    """
    Write the lock file pinning the enabled pipelines' catalogue datasets to
    ``datapkg.release``: only the items named by layers.csv ``datapkg_item``
    (whole datasets where a layer is resolved by filter). Commit the result;
    delete it to re-lock.
    """
    output:
        lock=DATAPKG["lock"],
    params:
        **_DATAPKG_LOCK_PARAMS,
    shell:
        """
        test -n "{params.base}" -a -n "{params.release}" -a -n "{params.selection}" || \
            (echo "datapkg_lock: set datapkg.base, datapkg.release and datapkg.enabled" && exit 1)
        irv-datapkg lock --base '{params.base}' --release '{params.release}' \
            {params.private_base} --output {output.lock} {params.selection}
        """


rule datapkg_lock_check:
    """
    CI: exit non-zero unless the lock file matches ``datapkg.release`` and
    locks every item the enabled pipelines need with its current checksum
    (`irv-datapkg lock ... --check`; writes nothing).
    """
    params:
        **_DATAPKG_LOCK_PARAMS,
        lock=DATAPKG["lock"],
    shell:
        """
        test -n "{params.base}" -a -n "{params.release}" -a -n "{params.selection}" || \
            (echo "datapkg_lock_check: set datapkg.base, datapkg.release and datapkg.enabled" && exit 1)
        irv-datapkg lock --check --base '{params.base}' --release '{params.release}' \
            {params.private_base} --output '{params.lock}' {params.selection}
        """
