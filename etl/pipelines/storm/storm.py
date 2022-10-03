"""
Helpers and Handlers for STORM (Tropical Cyclones)

For use within the Snakemake workflow
"""
import os
import sys
from typing import List

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helpers import get_logger, tiff_appears_valid


LOG = get_logger()


class HazardStorm:
    """
    Structure assumed:
        /top_storm_dir
            /present/.tiffs
            /future/.tiffs
    """

    def __init__(self):
        pass

    def get_filepaths(self, top_dir: str, limit=None) -> List[str]:
        """
        Collect and return a list of tiff filepaths to be processed

        ::arg top_dir str Top-level directory containing 'present' amd 'future' folders
        ::kwarg int limit Arbitrarily limit the number of collected files
        """
        file_paths = []
        if not os.path.exists(os.path.join(top_dir, "present")):
            LOG.error("missing directory `present` under %s", top_dir)
        if not os.path.exists(os.path.join(top_dir, "future")):
            LOG.error("missing directory `future` under %s", top_dir)
        for dirpath, _, filenames in os.walk(os.path.join(top_dir, "present")):
            for _file in filenames:
                fpath = os.path.join(dirpath, _file)
                if tiff_appears_valid(fpath) is True:
                    file_paths.append(fpath)
                else:
                    LOG.warning("tiff %s appears to be invalid", _file)

        for dirpath, _, filenames in os.walk(os.path.join(top_dir, "future")):
            for _file in filenames:
                fpath = os.path.join(dirpath, _file)
                if tiff_appears_valid(fpath) is True:
                    file_paths.append(fpath)
                else:
                    LOG.warning("tiff %s appears to be invalid", _file)
        LOG.debug("generated file_paths: %s", file_paths)
        return file_paths

    def file_metadata(self, filepath: str) -> dict:
        """
        Generate metadata for the given filepath, using the filename dictionary
            .../future/STORM_FIXED_RETURN_PERIODS_CMCC-CM2-VHR4_EP_4000_YR_RP.tif

        URL: https://data.4tu.nl/articles/dataset/STORM_climate_change_tropical_cyclone_wind_speed_return_periods/14510817/3

        The STORM_FIXED_RETURN_PERIOD datasets contain maximum wind speeds for a
            fixed set of return periods at 10 km resolution in every basin
            and for every climate model used here (see below).

        The STORM_FIXED_WIND_SPEED dataset contains return periods
            for a fixed set of maximum wind speeds at 10 km resolution in every ocean basin
        """
        meta = {
            "filename": "",
            "filepath": "",
            "epoch": "",
            "gcm": "",
            "rp": "",
            "wind_speed": "",
        }

        # File sub-dir process for epoch
        if "present" in filepath:
            meta["epoch"] = "2015"
        elif "future in filepath":
            meta["epoch"] = "2050"
        else:
            raise Exception("file does not conform (epoch): %s", filepath)
        # Parse fname
        filename = os.path.basename(filepath.lower())
        meta["filename"] = filename
        meta["filepath"] = filepath
        parsed = filename.split(".")[0].split("_")
        LOG.debug("parsed filename: %s", parsed)
        if parsed[1:4] == ["fixed", "wind", "speeds"]:
            meta["wind_speed"] = parsed[6]
            meta["rp"] = "none"
        elif parsed[1:4] == ["fixed", "return", "periods"]:
            meta["wind_speed"] = "none"
            meta["rp"] = parsed[6]
        else:
            raise Exception(
                "file does not conform (windspeed/rp was %s): %s", parsed[1:4], filepath
            )
            pass
        meta["gcm"] = parsed[4]
        LOG.debug("final meta: %s", meta)
        return meta

    # parse fname -> Dict / List
    # Ingest to API - Metadata
    # Ingest to terracotta sqlite DB (Python API)


if __name__ == "__main__":
    h = HazardStorm()
    file_paths = h.get_filepaths(
        "/home/dusted/code/oxford/infra-risk-vis/etl/raw_data/processed_data/input/storm"
    )
    for f in file_paths:
        h.file_metadata(f)
