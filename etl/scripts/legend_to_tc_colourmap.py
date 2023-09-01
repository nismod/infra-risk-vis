"""
Generate a Terracotta explicit_color_map for a given CSV file
"""

import csv
import argparse
import json

parser = argparse.ArgumentParser(description="Explicit ColorMap Generator")
parser.add_argument(
    "csv_fpath",
    type=str,
    help="Path to input (legend CSV file)",
)
parser.add_argument(
    "json_fpath",
    type=str,
    help="Path to output (legend JSON file)",
)


def main(
    csv_fpath,
    value_col="Value",
    red_col="R",
    green_col="G",
    blue_col="B",
    alpha=255,
    include_names=True,
    name_col="Label",
):
    with open(csv_fpath, mode="r", encoding="utf-8-sig") as csvfile:
        reader = csv.DictReader(csvfile)
        output = {}
        for line in reader:
            if include_names is True:
                value = {
                    "name": line[name_col],
                    "color": (
                        int(line[red_col]),
                        int(line[green_col]),
                        int(line[blue_col]),
                        alpha,
                    ),
                }
            else:
                value = (
                    int(line[red_col]),
                    int(line[green_col]),
                    int(line[blue_col]),
                    alpha,
                )
            output[int(line[value_col])] = value
    return output


if __name__ == "__main__":
    args = parser.parse_args()
    output = main(args.csv_fpath)

    with open(args.json_fpath, "w") as fptr:
        json.dump(output, fptr, indent=4)
