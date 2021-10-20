#!/bin/env python3

from argparse import ArgumentParser
import csv
import os
from pathlib import Path
import subprocess
import sys

this_directory = Path(__file__).parent.resolve()
vector_script_path = this_directory / 'prepare_vector.sh'


def run_single_processing(in_file_path: Path, out_file_path: Path, layer_name: str, output_layer_name: str, spatial_type: str, where_filter: str, **kwargs):
  print(f'Processing vector "{in_file_path}" -> "{out_file_path}"')
  command = f'{vector_script_path} "{in_file_path}" "{out_file_path}" "{output_layer_name}" "{spatial_type}" "{layer_name}" "{where_filter}"'
  print(f"Running command: {command}", flush=True)
  subprocess.run(command, shell=True, stdout=sys.stdout, stderr=sys.stderr)


def process_vector_datasets(raw: Path, out: Path):
  infrastructure_dir = raw / 'networks'
  csv_path = infrastructure_dir / 'network_layers.csv'
  assert csv_path.is_file(), f"{csv_path} is not a file"

  with csv_path.open() as f:
    reader = csv.DictReader(f)
    assert 'path' in reader.fieldnames
    assert 'layer_name' in reader.fieldnames
    assert 'spatial_type' in reader.fieldnames
    assert 'where_filter' in reader.fieldnames
    assert 'output_layer_name' in reader.fieldnames

    for row in reader:
      in_file_path = raw / row['path']
      output_layer_name = row['output_layer_name']
      out_file_path = out / f"{output_layer_name}.mbtiles"
      if os.path.exists(out_file_path) and (os.path.getmtime(in_file_path) < os.path.getmtime(out_file_path)):
        print("Skipping", out_file_path)
        continue
      run_single_processing(in_file_path, out_file_path, **row)


if __name__ == '__main__':
  parser = ArgumentParser(description='Converts all vector datasets to GeoJSON and then to MBTILES')
  parser.add_argument('--raw', type=Path, help='Root of the raw data directory. Assumes a file network_layers.csv exists in the dir.', required=True)
  parser.add_argument('--out', type=Path, help='Directory in which to store results of the processing', required=True)

  args = parser.parse_args()

  process_vector_datasets(args.raw.expanduser().resolve(), args.out.expanduser().resolve())
