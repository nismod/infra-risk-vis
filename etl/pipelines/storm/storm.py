"""
Helpers and Handlers for STORM (Tropical Cyclones)

For use within the Snakemake workflow
"""
import argparse
import csv
import os
import sys
from typing import List
from dataclasses import dataclass, asdict

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helpers import (
    count_hazard_csv_rows,
    file_key_from_fname,
    fname_from_fpath,
    get_logger,
    hazard_csv_exists,
    hazard_csv_valid,
    tiff_appears_valid,
)


LOG = get_logger()

parser = argparse.ArgumentParser(description="STORM Processor")
parser.add_argument(
    "source_dir",
    type=str,
    help="Source Directory for files",
    default=os.path.dirname(os.path.abspath(__file__)),
)
parser.add_argument(
    "--hazard_csv_fpath", type=str, help="path to hazard file", default=None
)


@dataclass
class StormMeta:
    """Metadata associated with an input file"""

    hazard: str
    key: str
    fname: str
    gcm: str
    rp: int


class StormFile:
    def __init__(self, fpath: str):
        self.source = "STORM"
        self.fpath = fpath
        self.fsize = None
        self.fname = fname_from_fpath(self.fpath)
        self.meta = self._parse_fname(self.fname)

    def _parse_fname(self, fname: str) -> StormMeta:
        """
        Parse the filepath for hazard.csv keys

        ::param fname st e.g. lange2020_hwmid-humidex_gfdl-esm2m_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005.nc4
            {MAINTAINER}_{MODEL}_{CLIMATE_FORCING}_{BIAS_ADJ}_{CLIMATE_SCENARIO(RCP)}_{SCO_ECO_SCENARIO}_{SENS_SCENARIO}_{VARIABLE}_{REGION}_{TIME_STEP}_{YEAR_START}_{YEAR_END}.nc4

        ::returns meta ISIMPExtremeHeatMeta
        """
        _key = file_key_from_fname(fname)
        parts = _key.split("_")
        # NOTE: hazard name matching front-end
        return StormMeta(
            hazard="cyclone", fname=fname, key=_key, gcm=parts[4], rp=int(parts[5])
        )


class HazardStorm:
    def __init__(self):
        self.hazard_csv_fieldnames = [
            "hazard",
            "path",
            "rp",
            "gcm",
            "key",
        ]

    def get_metadata(self, top_dir: str) -> List[StormFile]:
        """
        Collect and return a list of tiff filepaths to be processed

        ::arg top_dir str Top-level directory
        """
        input_files = []
        for dirpath, _, filenames in os.walk(top_dir):
            for _file in filenames:
                if "STORM_FIXED_RETURN_PERIODS" in _file:
                    fpath = os.path.join(dirpath, _file)
                    if tiff_appears_valid(fpath) is True:
                        input_files.append(self.file_metadata(fpath))
                    else:
                        LOG.warning("tiff %s appears to be invalid", _file)
        LOG.debug("generated file data for %s input_files", len(input_files))
        return input_files

    def file_metadata(self, fpath: str) -> StormFile:
        """
        Generate metadata for the given fpath, using the filename
            STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_EP_4000_YR_RP.tif
        NOTE: Assumes we are working with post-processed STORM data (mosaiced to Global)
        """
        storm_file = StormFile(fpath)
        LOG.debug("final storm meta: %s", asdict(storm_file.meta))
        return storm_file

    def append_hazard_csv(
        self, input_files: List[StormFile], hazard_csv_fpath: str
    ) -> None:
        """
        Append the meta from parsed files to the hazard csv
        """
        # Generate file if req
        if not hazard_csv_exists(hazard_csv_fpath):
            with open(hazard_csv_fpath, "w") as csvfile:
                LOG.warning(" No hazard csv found writing new at: %s", hazard_csv_fpath)
                # Generate a new file
                writer = csv.DictWriter(csvfile, fieldnames=self.hazard_csv_fieldnames)
                writer.writeheader()
        # Check CSV is valid (headers etc)
        _ = hazard_csv_valid(hazard_csv_fpath, self.hazard_csv_fieldnames)
        # Dump meta
        count_before = count_hazard_csv_rows(hazard_csv_fpath)
        with open(hazard_csv_fpath, "a") as csvfile:
            # Setup the writer
            writer = csv.DictWriter(csvfile, fieldnames=self.hazard_csv_fieldnames)
            # Dump meta
            count_outputs = 0
            for input_file in input_files:
                row = {}
                row["hazard"] = input_file.meta.hazard
                row["rp"] = input_file.meta.rp
                row["gcm"] = input_file.meta.gcm
                row["key"] = input_file.meta.key
                row["path"] = input_file.fpath
                writer.writerow(row)
                count_outputs += 1
        count_after = count_hazard_csv_rows(hazard_csv_fpath)
        LOG.info(
            "Count before: %s, Count after: %s, number input files processed: %s",
            count_before - 1,
            count_after,
            len(input_files),
        )


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    processor = HazardStorm()
    input_files = processor.get_metadata(args.source_dir)
    # Generate the CSV
    processor.append_hazard_csv(input_files, args.hazard_csv_fpath)
    LOG.info("Done.")
