"""Fetch one pipeline layer from the irv-datapkg catalogue to raster/raw/.

Called by the ``fetch_datapkg`` rule (rules/datapkg.smk); can be run by hand:

    python scripts/datapkg_fetch.py --pipeline aqueduct \
        --key inunriver_rcp4p5_00000NorESM1-M_2030_rp00100 \
        --output raster/raw/aqueduct/inunriver_rcp4p5_00000NorESM1-M_2030_rp00100.tif \
        --base s3://bucket/irv-datapkg --release 2026.1
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import datapkg_mapping as dm  # noqa: E402


def fetch_layer(pipeline, key, output, base=None, release=None, lock=None, private_base=None,
                cache_dir=None, spec=None, root="."):
    """Resolve ``pipeline``/``key`` to a catalogue item and fetch it to ``output``,
    verifying sha256 (against the lock file when one is used)."""
    if spec is None:
        mapping = dm.load_mapping(pipeline, root)
        spec = dm.item_spec(dm.layer_row(dm.read_layers(pipeline, root), key), mapping)
    catalog = dm.open_catalog(base, release, lock, private_base, cache_dir)
    item = dm.resolve(catalog, spec)
    if spec.get("sha256") and item.asset.sha256 != spec["sha256"]:
        raise dm.MappingError(
            f"{item.dataset}/{item.id}: catalogue sha256 differs from the lock file"
        )
    output = Path(output)
    tmp = output.with_name(output.name + ".part")
    catalog.fetch(item, tmp)
    tmp.replace(output)
    return item, output


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
    item, path = fetch_layer(
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
    print(f"{item.dataset}/{item.version}/{item.id} -> {path}")


if __name__ == "__main__":
    main()
