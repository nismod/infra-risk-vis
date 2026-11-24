"""Fetch one pipeline layer from the irv-datapkg catalogue to raster/raw/.

Called by the ``fetch_datapkg`` rule (rules/datapkg.smk); can be run by hand:

    python scripts/datapkg_fetch.py --pipeline aqueduct \
        --key inunriver_rcp4p5_00000NorESM1-M_2030_rp00100 \
        --output raster/raw/aqueduct/inunriver_rcp4p5_00000NorESM1-M_2030_rp00100.tif \
        --base s3://bucket/irv-datapkg --release 2026.1

With ``--lock`` (an existing lock file), ``--base`` overrides the lock's
recorded base, e.g. to read from a mirror; ``--release`` is ignored.
"""

import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import datapkg_mapping as dm  # noqa: E402


def fetch_layer(pipeline, key, output, base=None, release=None, lock=None, private_base=None,
                cache_dir=None, spec=None, root="."):
    """Resolve ``pipeline``/``key`` to a catalogue asset and fetch it to ``output``.

    The client verifies size and sha256 while streaming and renames into
    place atomically. With a lock file and a spec naming an item, the file is
    fetched from the lock entry alone (href, sha256, size: no STAC reads);
    otherwise the item is resolved in the catalogue, and from a lock the
    client refuses an item whose checksum differs from it.

    Without ``cache_dir`` (or ``$IRV_DATAPKG_CACHE``) the file is streamed
    straight to ``output`` (no second copy on disk); with one, it is cached
    there and hardlinked to ``output`` (read-only, shared with the cache).

    Returns ``(asset, output)``.
    """
    if spec is None:
        mapping = dm.load_mapping(pipeline, root)
        spec = dm.item_spec(dm.layer_row(dm.read_layers(pipeline, root), key), mapping)
    cache_dir = cache_dir or os.environ.get("IRV_DATAPKG_CACHE")
    catalog = dm.open_catalog(base, release, lock, private_base, cache_dir)
    how = dict(link="hardlink", cache=True) if cache_dir else dict(cache=False)
    output = Path(output)
    locked = catalog.lockfile
    if locked is not None and "item" in spec and spec["dataset"] in locked.datasets:
        asset = catalog.locked_asset(spec["dataset"], spec["item"])
        catalog.fetch_locked(spec["dataset"], spec["item"], output, **how)
        return asset, output
    item = dm.resolve(catalog, spec)
    if locked is None and spec.get("sha256") and item.asset.sha256 != spec["sha256"]:
        # only reachable by hand (the rule adds sha256 from a lock file)
        raise dm.MappingError(f"{item.dataset}/{item.id}: catalogue sha256 differs from the spec")
    catalog.fetch(item, output, **how)
    return item.asset, output


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pipeline", required=True)
    p.add_argument("--key", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--spec", help="JSON item spec (default: resolved from layers.csv)")
    p.add_argument("--base")
    p.add_argument("--private-base")
    p.add_argument("--release")
    p.add_argument("--lock")
    p.add_argument("--cache-dir")
    args = p.parse_args(argv)
    asset, path = fetch_layer(
        args.pipeline,
        args.key,
        args.output,
        base=args.base,
        release=args.release,
        lock=args.lock,
        private_base=args.private_base,
        cache_dir=args.cache_dir,
        spec=json.loads(args.spec) if args.spec else None,
    )
    print(f"{asset.dataset}/{asset.version}/{asset.item_id} -> {path}")


if __name__ == "__main__":
    main()
