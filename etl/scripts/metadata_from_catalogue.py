"""Fill a pipeline's metadata.json ``license`` and ``description`` from its
irv-datapkg STAC collection.

``domain``, ``name``, ``group`` and ``keys`` (and colour maps, which live in
separate files) stay local. ``description`` takes the collection's
``sci:citation`` (+ DOI) by default, because raster_tile_sources has no
citation column and the frontend shows ``description`` as the source credit;
``--description description`` uses the collection description instead.

    python scripts/metadata_from_catalogue.py aqueduct --base s3://.../irv-datapkg
    python scripts/metadata_from_catalogue.py aqueduct --check   # CI: exit 1 if stale

Catalogue options default to the ``datapkg`` section of config.yml.
"""

import argparse
import json
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).parent))

import datapkg_mapping as dm  # noqa: E402


def update_metadata(pipeline, catalog, root=".", description="citation"):
    mapping = dm.load_mapping(pipeline, root)
    path = dm.pipeline_dir(pipeline, root) / "metadata.json"
    metadata = json.loads(path.read_text())
    collection = catalog.collection(mapping["dataset"])
    return path, dm.metadata_from_collection(metadata, collection, description)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pipeline")
    p.add_argument("--root", default=".")
    p.add_argument("--config", default=None)
    p.add_argument("--base")
    p.add_argument("--private-base")
    p.add_argument("--release")
    p.add_argument("--lock")
    p.add_argument("--description", choices=["citation", "description"], default="citation")
    p.add_argument("--check", action="store_true")
    args = p.parse_args(argv)

    config_path = Path(args.config or Path(args.root, "config.yml"))
    config = yaml.safe_load(config_path.read_text()) if config_path.exists() else {}
    cfg = dm.datapkg_config(config)
    lock = args.lock or (str(Path(args.root, cfg["lock"])) if cfg["lock"] else None)
    catalog = dm.open_catalog(
        base=args.base or cfg["base"],
        release=args.release or cfg["release"],
        lock=lock if (args.lock or not args.base) else None,
        private_base=args.private_base or cfg["private_base"],
        cache_dir=cfg["cache_dir"],
    )
    path, metadata = update_metadata(args.pipeline, catalog, args.root, args.description)
    text = dm.dump_metadata(metadata)
    if args.check:
        if path.read_text() != text:
            print(f"{path} differs from the catalogue collection", file=sys.stderr)
            return 1
        return 0
    path.write_text(text)
    print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
