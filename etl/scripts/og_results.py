"""Read open-gira results published as an irv-datapkg results run.

open-gira publishes the outputs infra-risk-vis displays to
``{BASE}/results/open-gira/{run}/`` (irv-datapkg ``docs/catalogue-spec.md``
section 8): a static STAC catalogue with one collection per layer
(``road_edges_damages``, ``rail_edges``, ``rail_nodes``, ``nbs_opportunities``,
``admin_boundaries``, ``hydrobasins_gadm``, ...) whose only version is the run
id, holding GeoParquet 1.1 files. The vector pipelines (osm_roads, osm_rail,
nbs-adaptation) use this module to fetch those layers with the irv-datapkg
client instead of reading copies under ``/opt/etl/raw_data/`` or
``/mnt/linux-filestore/...``.

Configured per pipeline by an ``open_gira_results`` block::

    open_gira_results:
      base: s3://nismod/irv-datapkg      # irv-datapkg {BASE}
      run: 20261124-fc04e93-r1           # published run id
      lock: null                         # optional results lock (pins checksums)
      cache_dir: null                    # download cache (default $IRV_DATAPKG_CACHE)
      storage_options: {}

With neither ``run`` nor an existing ``lock`` the pipelines keep their legacy
paths. Command line (run from etl/)::

    python scripts/og_results.py fetch --base B --run R --layer rail_edges \\
        --filter osm_extract=planet-latest --output raw_data/x.geoparquet
    python scripts/og_results.py lock --base B --run R --layer rail_edges \\
        --output pipelines/osm_rail/og-results.lock.yaml
"""

import argparse
import os
import posixpath
import sys
from pathlib import Path

import yaml

PRODUCER = "open-gira"
DEFAULTS = {
    "base": None,
    "run": None,
    "lock": None,
    "cache_dir": None,
    "storage_options": {},
    # fetched files go to {dest}/{run}/...
    "dest": "raw_data/open-gira",
}


class ResultsError(Exception):
    pass


def _opt(value):
    """``--config`` passes nested nulls as strings."""
    if value is None or (isinstance(value, str) and value.strip() in ("", "null", "None", "~")):
        return None
    return value


def results_config(config: dict) -> dict:
    """The ``open_gira_results`` block with defaults filled in."""
    cfg = {**DEFAULTS, **((config or {}).get("open_gira_results") or {})}
    for key in ("base", "run", "lock", "cache_dir"):
        cfg[key] = _opt(cfg[key])
    if cfg["lock"] and os.path.exists(cfg["lock"]) and not cfg["run"]:
        cfg["run"] = (yaml.safe_load(Path(cfg["lock"]).read_text()) or {}).get("run")
    return cfg


def active(cfg: dict) -> bool:
    """True if the pipeline should read from a published run."""
    return bool(cfg["run"] or (cfg["lock"] and os.path.exists(cfg["lock"])))


def run_root(base: str, run: str) -> str:
    """``{BASE}/results/open-gira/{run}``."""
    return str(base).rstrip("/") + f"/results/{PRODUCER}/{run}"


def parse_filters(text) -> dict:
    """``"osm_extract=planet-latest;hazard=aqueduct-river"`` -> dict. Digits
    stay strings: the client compares numbers numerically either way."""
    if text is None or (isinstance(text, float) and text != text) or not str(text).strip():
        return {}
    if isinstance(text, dict):
        return dict(text)
    out = {}
    for part in str(text).split(";"):
        if not part.strip():
            continue
        key, sep, value = part.partition("=")
        if not sep:
            raise ResultsError(f"filter {part!r} is not key=value")
        out[key.strip()] = value.strip()
    return out


def open_run(cfg: dict):
    """The run's catalogue: from the lock file if it exists, else from
    ``base`` + ``run``. Checks the run root is an open-gira results run."""
    from irv_datapkg.catalog import Catalog
    from irv_datapkg.catalog import _fs

    kwargs = {"storage_options": cfg.get("storage_options") or {}}
    if cfg.get("cache_dir") or os.environ.get("IRV_DATAPKG_CACHE"):
        kwargs["cache_dir"] = cfg.get("cache_dir") or os.environ["IRV_DATAPKG_CACHE"]
    lock = cfg.get("lock")
    if lock and os.path.exists(lock):
        base = run_root(cfg["base"], cfg["run"]) if cfg.get("base") and cfg.get("run") else None
        catalog = Catalog.from_lock(lock, base=base, **kwargs)
        run = catalog.lockfile.data.get("run")
        if not run:
            raise ResultsError(f"{lock} is not a results lock (no 'run')")
        if cfg.get("run") and cfg["run"] != run:
            raise ResultsError(f"{lock} pins run {run}, config asks for {cfg['run']}")
        return catalog
    if not (cfg.get("base") and cfg.get("run")):
        raise ResultsError("open_gira_results: set base and run (or an existing lock)")
    root = run_root(cfg["base"], cfg["run"])
    catalog = Catalog.open(root, **kwargs)
    try:
        doc = _fs.read_json(_fs.join(catalog.base, "catalog.json"), kwargs["storage_options"])
    except FileNotFoundError:
        raise ResultsError(f"no published run at {root}") from None
    if doc.get("irv:producer") != PRODUCER or doc.get("irv:run") != cfg["run"]:
        raise ResultsError(
            f"{root}/catalog.json is not open-gira run {cfg['run']} "
            f"(irv:producer={doc.get('irv:producer')!r}, irv:run={doc.get('irv:run')!r})"
        )
    return catalog


def _place(catalog, item, dest: Path, cache: bool) -> Path:
    """Fetch one item's data asset to ``dest`` (sha256-verified). With a cache,
    it is downloaded once and hardlinked; otherwise streamed straight to dest."""
    how = {"link": "hardlink", "cache": True} if cache else {"cache": False}
    if catalog.lockfile is not None and item.dataset in catalog.lockfile.datasets:
        return catalog.fetch_locked(item.dataset, item.id, dest, **how)
    return catalog.fetch(item, dest, **how)


def use_cache(cfg: dict) -> bool:
    return bool(cfg.get("cache_dir") or os.environ.get("IRV_DATAPKG_CACHE"))


def select(catalog, layer: str, filters: dict) -> list:
    items = catalog.items(layer, **filters)
    if not items:
        raise ResultsError(f"{layer}: no items match {filters or 'any'}")
    return sorted(items, key=lambda i: i.id)


def fetch_layer(catalog, layer: str, output, filters=None, single=False, cache=False) -> list:
    """Fetch a layer's items (matching ``filters``) and return their paths.

    ``single``: exactly one item must match; it is written to ``output``
    (a file). Otherwise ``output`` is a directory holding
    ``{item id}.geoparquet`` per item (the osm_roads loader walks it)."""
    items = select(catalog, layer, parse_filters(filters))
    output = Path(output)
    if single:
        if len(items) != 1:
            raise ResultsError(
                f"{layer}: {len(items)} items match {filters or 'any'} "
                f"({', '.join(i.id for i in items)}); add filters to pick one"
            )
        return [_place(catalog, items[0], output, cache)]
    output.mkdir(parents=True, exist_ok=True)
    return [_place(catalog, item, output / f"{item.id}.geoparquet", cache) for item in items]


def lock(cfg: dict, layers: dict, path) -> Path:
    """Write a results lock (catalogue-spec section 8) for ``layers``
    ({layer: filters}): the run, and each selected item's href and sha256."""
    catalog = open_run({**cfg, "lock": None})
    datasets = {}
    for layer, filters in layers.items():
        entry = datasets.setdefault(
            layer,
            {"version": catalog.version(layer), "access": catalog.access(layer), "items": {}},
        )
        dims = catalog.dimensions(layer)
        for item in select(catalog, layer, parse_filters(filters)):
            entry["items"][item.id] = {
                key: {
                    "href": catalog.relative_href(a),
                    "sha256": a.sha256,
                    "size": a.size,
                    **({"source_filename": a.source_filename} if a.source_filename else {}),
                }
                for key, a in item.assets.items()
            }
            if dims:
                entry.setdefault("properties", {})[item.id] = {
                    d: item.properties.get(d) for d in dims
                }
    data = {
        "lock_version": 1,
        "base": run_root(cfg["base"], cfg["run"]),
        "run": cfg["run"],
        "datasets": datasets,
    }
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(data, sort_keys=False))
    return path


def fetched_path(cfg: dict, ref: str, single: bool) -> str:
    """Where a pipeline fetches a network layer ``ref`` to."""
    return posixpath.join(cfg["dest"], str(cfg["run"]), ref + (".geoparquet" if single else ""))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("fetch", "lock"):
        s = sub.add_parser(name)
        s.add_argument("--base")
        s.add_argument("--run")
        s.add_argument("--layer", required=True, action="append")
        s.add_argument("--filter", action="append", default=[],
                       help="key=value[;key=value], applied to every --layer")
        s.add_argument("--output", required=True)
        if name == "fetch":
            s.add_argument("--lock")
            s.add_argument("--single", action="store_true",
                           help="exactly one item: write it to --output (a file)")
        s.add_argument("--cache-dir")
    args = p.parse_args(argv)
    filters = {}
    for f in args.filter:
        filters.update(parse_filters(f))
    cfg = results_config({"open_gira_results": {
        "base": args.base, "run": args.run, "lock": getattr(args, "lock", None),
        "cache_dir": args.cache_dir}})
    if args.command == "lock":
        print(lock(cfg, {layer: filters for layer in args.layer}, args.output))
        return 0
    if len(args.layer) != 1:
        p.error("fetch takes one --layer")
    for path in fetch_layer(open_run(cfg), args.layer[0], args.output, filters, args.single,
                            cache=use_cache(cfg)):
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())


__all__ = ["active", "fetch_layer", "lock", "open_run", "results_config", "run_root"]
