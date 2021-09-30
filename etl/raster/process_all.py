#!/bin/env python3

from argparse import ArgumentParser
import csv
from pathlib import Path
import subprocess

this_directory = Path(__file__).parent.resolve()
raster_script_path = this_directory / 'prepare_raster.sh'

def run_single_processing(in_file_path: Path, out_file_path: Path):
  print(f'Processing raster "{in_file_path}" -> "{out_file_path}"')
  subprocess.call(f'{raster_script_path} "{in_file_path}" "{out_file_path}"', shell=True)


def process_raster_datasets(raw: Path, out: Path, type: str = None):
  csv_path = raw / 'hazard_layers.csv'
  assert csv_path.is_file(), f"{csv_path} is not a file"

  with csv_path.open() as f:
    reader = csv.DictReader(f)
    assert 'path' in reader.fieldnames
    assert 'key' in reader.fieldnames
    if type is not None:
      assert 'hazard' in reader.fieldnames
    for row in reader:
      if type is None or row['hazard'] == type:
        in_file_path = raw / row['path']
        file_key = row['key'].replace('.', 'x') # replace dots with 'x' because dots don't work with terracotta
        out_file_path = out / f"{file_key}.tif"
        run_single_processing(in_file_path, out_file_path)


if __name__ == '__main__':
  parser = ArgumentParser(description='Converts all raster datasets to Cloud-Optimized GeoTIFFs')
  parser.add_argument('--raw', type=Path, help='Root of the raw data directory. Assumes a file hazard_layers.csv exists in the dir.', required=True)
  parser.add_argument('--out', type=Path, help='Directory in which to store results of the processing', required=True)
  parser.add_argument('--type', type=str, help='Type of hazard layer to process. Will ignore other types of layers.', default=None)

  args = parser.parse_args()

  process_raster_datasets(args.raw.expanduser().resolve(), args.out.expanduser().resolve(), args.type)

