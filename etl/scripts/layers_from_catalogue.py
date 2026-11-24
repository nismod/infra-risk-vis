"""Generate or refresh a pipeline's layers.csv from irv-datapkg STAC items.

Refresh (default): match every existing row of pipelines/{pipeline}/layers.csv
to exactly one catalogue item using pipelines/{pipeline}/datapkg.yml, and
record its id in the ``datapkg_item`` column. Filenames (the raster keys used
by terracotta and the Zarr stack) and all other columns are kept.

    python scripts/layers_from_catalogue.py aqueduct --base s3://.../irv-datapkg --release 2026.1

``--add-missing`` also appends a row for every item no row matched, deriving
key columns from item properties (reverse of the datapkg.yml value maps) and
``filename``/``columns`` from its templates. Run it on a pipeline with no
layers.csv to generate one from scratch. ``--check`` writes nothing and exits
non-zero if the file would change or any row is unmatched/ambiguous.

Catalogue options default to the ``datapkg`` section of config.yml.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent))

import datapkg_mapping as dm  # noqa: E402


def generate(pipeline, catalog, root=".", add_missing=False):
    mapping = dm.load_mapping(pipeline, root)
    path = dm.pipeline_dir(pipeline, root) / "layers.csv"
    if path.exists():
        layers = dm.read_layers(pipeline, root)
    else:
        columns = ["filename", *mapping["properties"], *mapping["columns"]]
        layers = pd.DataFrame(columns=list(dict.fromkeys(columns)))
        add_missing = True
    items = catalog.items(mapping["dataset"])
    return dm.refresh_layers(layers, items, mapping, add_missing=add_missing)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("pipeline")
    p.add_argument("--root", default=".", help="etl directory (default: .)")
    p.add_argument("--config", default=None, help="config.yml (default: {root}/config.yml)")
    p.add_argument("--base")
    p.add_argument("--private-base")
    p.add_argument("--release")
    p.add_argument("--lock")
    p.add_argument("--add-missing", action="store_true")
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

    layers, report = generate(args.pipeline, catalog, args.root, args.add_missing)
    for name in report["unmatched"]:
        print(f"unmatched row: {name}", file=sys.stderr)
    for name, ids in report["ambiguous"]:
        print(f"ambiguous row: {name} -> {', '.join(ids)}", file=sys.stderr)
    for item_id in report["unused"]:
        print(f"item not in layers.csv: {item_id}", file=sys.stderr)
    for name in report["new"]:
        print(f"new row: {name}", file=sys.stderr)

    path = dm.pipeline_dir(args.pipeline, args.root) / "layers.csv"
    text = layers.to_csv(index=False, lineterminator="\n")
    failed = bool(report["unmatched"] or report["ambiguous"])
    if args.check:
        current = path.read_text() if path.exists() else ""
        if current != text:
            print(f"{path} is out of date", file=sys.stderr)
            failed = True
    else:
        path.write_text(text)
        print(path)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
