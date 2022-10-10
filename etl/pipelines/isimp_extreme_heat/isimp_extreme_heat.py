"""
Data Set Downloaders
"""

import os
import sys
import csv
from typing import List
import logging
import argparse

import numpy as np
from netCDF4 import Dataset
from osgeo import gdal
from py import process

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from pipelines.helpers import _download_file, file_key_from_fname

logging.basicConfig(format="%(asctime)s %(levelname)s:%(message)s", level=logging.DEBUG)


parser = argparse.ArgumentParser(description="ISIMP Extreme Heat Download / Processor")
parser.add_argument(
    "download_dir",
    type=str,
    help="Download directory for image files, defaults to script directory",
    default=os.path.dirname(os.path.abspath(__file__)),
)
parser.add_argument(
    "--output_occurrence_data_directory",
    type=str,
    help="Directory in which to put processed heatwave occurrence TIFF files",
    default=os.path.dirname(os.path.abspath(__file__)),
)
parser.add_argument(
    "--output_exposure_data_directory",
    type=str,
    help="Directory in which to put processed heatwave exposure TIFF files",
    default=os.path.dirname(os.path.abspath(__file__)),
)
parser.add_argument(
    "--hazard_csv_fpath", type=str, help="path to cumulative hazard file", default=None
)
parser.add_argument(
    "--limit_files",
    type=int,
    default=None,
    help="Limit the number of files processed arbitrarily (for testing) - True / False",
)
parser.add_argument(
    "--log_meta",
    type=str,
    default="False",
    help="Logout metadata associated with files",
)


class HazardISIMPExtremeHeat:
    """
    ISIMP Extreme Heat
    """

    def __init__(
        self,
        list_file: str = "https://github.com/nismod/open-gira/files/9488963/GRII_open_hazard_links.csv",
        download_dir: str = None,
        hazard_csv_fpath: str = None,
        world_pop_2020_05deg_fpath: str = "worldpop_2020_resamp_0x5deg.tif",
    ):
        self.list_file = list_file
        self.hazard_csv_fpath = hazard_csv_fpath
        self.download_dir = download_dir
        self.world_pop_2020_05deg_fpath = world_pop_2020_05deg_fpath
        self.hazard_name = "extreme_heat"
        self.tmp_dir = "/tmp"  # For list file processing
        self.hazard_csv_fieldnames = (
            [
                "hazard",
                "path",
                "rcp",
                "epoch",
                "gcm",
                "key",
            ],
        )
        self.epoch_bins = {
            "baseline": {
                "bin_start": 1966,
                "bin_end": 2005,
                "file_year_start": 1861,
                "file_year_end": 2005,
            },
            "2030": {
                "bin_start": 2010,
                "bin_end": 2049,
                "file_year_start": 2006,
                "file_year_end": 2099,
            },
            "2050": {
                "bin_start": 2030,
                "bin_end": 2069,
                "file_year_start": 2006,
                "file_year_end": 2099,
            },
            "2080": {
                "bin_start": 2060,
                "bin_end": 2099,
                "file_year_start": 2006,
                "file_year_end": 2099,
            },
        }
        self.climate_scenarios_to_process = ["historical", "rcp26", "rcp60"]
        self.time_spans_to_process = [[1861, 2005], [2006, 2099]]
        if not os.path.exists(self.world_pop_2020_05deg_fpath):
            logging.warning(
                "world pop data does not exist at given path - exposure calcs will fail: %s",
                self.world_pop_2020_05deg_fpath,
            )

    @staticmethod
    def nc_fname_from_url(url: str) -> str:
        """
        Parse filename for nc file from url
        """
        return url.split("/")[-1]

    def parse_fname(self, fname: str) -> dict:
        """
        Parse the filename for hazard.csv keys

        ::param fname st e.g. lange2020_hwmid-humidex_gfdl-esm2m_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005.nc4
            {MAINTAINER}_{MODEL}_{CLIMATE_FORCING}_{BIAS_ADJ}_{CLIMATE_SCENARIO(RCP)}_{SCO_ECO_SCENARIO}_{SENS_SCENARIO}_{VARIABLE}_{REGION}_{TIME_STEP}_{YEAR_START}_{YEAR_END}.nc4

        ::returns meta dict
            {   'filename': 'lange2020_hwmid-humidex_gfdl-esm2m_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005.nc4',
                'meta': { ...parsed variables in name  }
            }
        """
        basename = fname.split(".")[0]
        parts = basename.split("_")
        return {
            "hazard": self.hazard_name,
            "key": basename,
            "maintainer": parts[0],
            "model": parts[1],
            "climate_forcing": parts[2],
            "bias_adj": parts[3],
            "climate_scenario": parts[4],
            "socio_econ_scenario": parts[5],
            "sensitivity_scenario": parts[6],
            "variable": parts[7],
            "region": parts[8],
            "time_step": parts[9],
            "year_start": int(parts[10]),
            "year_end": int(parts[11]),
        }

    def fetch_nc_paths(self) -> List[dict]:
        """
        Fetch the NC CSV and parse out the ISIMP Extreme Heat paths
        """
        # fetch list to download dir then remove
        local_list_file = os.path.join(self.tmp_dir, "extreme_heat_list.csv")
        try:
            fsize = _download_file(self.list_file, local_list_file)
            if not fsize:
                raise FileNotFoundError("failed to download list file")
            with open(local_list_file, "r") as csvfile:
                reader = csv.DictReader(csvfile)
                files_meta = list(reader)
                # # Remove extraneous types
                files_meta = [
                    item
                    for item in files_meta
                    if item["\ufeffhazard"] == self.hazard_name
                ]
                # Append the filenames and parse filename to meta
                for item in files_meta:
                    item["fname"] = self.nc_fname_from_url(item["url"])
                    item["path"] = ""
                    item["meta"] = self.parse_fname(item["fname"])
                # Remove climate scenarios and years we do not require
                files_meta = [
                    item
                    for item in files_meta
                    if item["meta"]["climate_scenario"]
                    in self.climate_scenarios_to_process
                    and [item["meta"]["year_start"], item["meta"]["year_end"]]
                    in self.time_spans_to_process
                ]
                return files_meta
        except:
            raise
        finally:
            os.remove(local_list_file)

    def download_files(self, files_meta: List[dict], limit_files=None) -> List[dict]:
        """
        Download the files in the given meta

        ::returns files_meta - updated with local file paths
        """
        for idx, file_meta in enumerate(files_meta):
            logging.debug("Downloading %s", file_meta["url"])
            try:
                fpath = os.path.join(self.download_dir, file_meta["fname"])
                # Skip if file already exists
                if os.path.exists(fpath):
                    logging.info(
                        "file already exists locally - skipping download: %s",
                        file_meta["fname"],
                    )
                    fsize = os.path.getsize(fpath)
                else:
                    fsize = _download_file(file_meta["url"], fpath)
                    if not fsize:
                        raise Exception(
                            "download resulted in invalid file: {}".format(
                                file_meta["url"]
                            )
                        )
                file_meta["fsize"] = fsize
                file_meta["path"] = fpath
                if limit_files and idx >= limit_files:
                    logging.info("download aborting due to limit_files=%s", limit_files)
                    break
            except Exception as err:
                logging.exception("")
        return files_meta

    def np_to_geotiff(
        self,
        data: np.ndarray,
        output_fpath: str,
        width=720,
        height=360,
        ul_x=-180.0,
        ul_y=90.0,
        pixel_size=0.5,
    ) -> None:
        """
        Dump np array to GeoTiff
        """
        drv = gdal.GetDriverByName("GTiff")
        ds = drv.Create(output_fpath, width, height, 1, gdal.GDT_Float32)
        ds.SetGeoTransform((ul_x, pixel_size, 0, ul_y, 0, -pixel_size))
        ds.GetRasterBand(1).WriteArray(data)

    def occurrence_tiff_fname_from_meta(self, file_meta: dict, epoch: str) -> str:
        """
        Generate a name for a processed tiff using the files metadata
        """
        basename = file_key_from_fname(file_meta["fname"])
        return f"{basename}_{epoch}_occurence.tif"

    def exposure_tiff_fname_from_occurrence(self, occurrence_tif_fname: str) -> str:
        """
        Generate a name for a exposure tif from occurrence
        """
        basename = file_key_from_fname(occurrence_tif_fname)
        fname = basename.replace("occurrence", "exposure")
        return f"{fname}.tif"

    def generate_popn_exposure_tiff(
        self, temporally_averaged_tiff_fpath: str, output_dir: str
    ) -> str:
        """
        Generate a TIFF file of world population x heatwave occurrence
        """
        logging.info(
            "generating popn exposure tiff - input file: %s",
            temporally_averaged_tiff_fpath,
        )
        # Load world pop
        raster = gdal.Open(self.world_pop_2020_05deg_fpath)
        worldpop_data = raster.ReadAsArray()

        # Load Temporally averaged TIFF
        raster = gdal.Open(temporally_averaged_tiff_fpath)
        heat_data = raster.ReadAsArray()

        # Sanity check the shape - should both be single band and of same shape
        assert (
            worldpop_data.shape == heat_data.shape
        ), f"worldpop and heatwave tiffs are of differing shape - aborting ({ worldpop_data.shape} vs {heat_data.shape})"

        # Do exposure calc
        exposure = heat_data * worldpop_data
        output_fname = self.exposure_tiff_fname_from_occurrence(
            os.path.basename(temporally_averaged_tiff_fpath)
        )
        output_fpath = os.path.join(output_dir, output_fname)
        # Save
        self.np_to_geotiff(exposure, output_fpath)
        return output_fpath

    def generate_temporally_averaged_tiff(
        self,
        file_meta: dict,
        output_epoch: str,
        output_dir: str,
        averaging_year_start: int,
        averaging_year_end: int,
    ) -> str:
        """
        Generate a TIFF file using temporal averaging for the given time window
        """
        logging.info(
            "generating temporally averaged tiff - input file: %s epoch: %s: averaging_year_start: %s, averaging_year_end: %s",
            file_meta["fname"],
            output_epoch,
            averaging_year_start,
            averaging_year_end,
        )
        # Load nc4
        dataset = Dataset(file_meta["path"], "r", format="NETCDF4")
        # Load data
        leh_np = dataset.variables["leh"][:]
        # Calculate temporal slice idx
        start_idx = averaging_year_start - file_meta["meta"]["year_start"]
        end_idx = averaging_year_end - file_meta["meta"]["year_start"]
        logging.debug(
            "slicing file %s at idx %s (yr %s) to %s (yr %s)",
            file_meta["fname"],
            start_idx,
            file_meta["meta"]["year_start"] + start_idx,
            end_idx,
            file_meta["meta"]["year_start"] + end_idx,
        )
        # Generate mean
        leh_np_mean = np.mean(leh_np[start_idx:end_idx, :, :], axis=(0))
        logging.debug(
            "resulting average shape: %s, min: %s, max: %s",
            leh_np_mean.shape,
            leh_np_mean.min(),
            leh_np_mean.max(),
        )
        # Write Occurrence GeoTiff
        output_fname = self.occurrence_tiff_fname_from_meta(file_meta, output_epoch)
        output_fpath = os.path.join(output_dir, output_fname)
        logging.debug("writing occurrence geotiff to: %s", output_fpath)
        self.np_to_geotiff(leh_np_mean, output_fpath)
        logging.debug(
            "wrote occurrence geotiff of size %s to: %s",
            os.path.getsize(output_fpath),
            output_fpath,
        )
        return output_fpath

    def nc4_to_average_tiffs(self, files_meta: dict, output_dir: str) -> List[dict]:
        """
        Temporal averaging for netCDF4 files contained in files meta
            Ref: Dr Mark bernhofen

        Note we only process the given climate senario files

        ::param files_meta List[dict] - metadata generated from list_url files
        ::param output_dir str - Processing output directory

        ::return result_filenames List[dict]
        """
        output_fpaths = []
        for file_meta in files_meta:
            # In case some eextraneous files are present:
            if (
                file_meta["meta"]["climate_scenario"]
                not in self.climate_scenarios_to_process
            ):
                continue
            if file_meta["meta"]["climate_scenario"] == "historical":
                # Sanity check file meta start_ends
                assert (
                    file_meta["meta"]["year_start"]
                    == self.epoch_bins["baseline"]["file_year_start"]
                )
                assert (
                    file_meta["meta"]["year_end"]
                    == self.epoch_bins["baseline"]["file_year_end"]
                )
                try:
                    # Special case for historical
                    tiff_fpath = self.generate_temporally_averaged_tiff(
                        file_meta,
                        "baseline",
                        output_dir,
                        self.epoch_bins["baseline"]["bin_start"],
                        self.epoch_bins["baseline"]["bin_end"],
                    )
                except Exception:
                    logging.exception("")
                output_fpaths.append(tiff_fpath)
            else:
                for epoch, bin_range in self.epoch_bins.items():
                    if epoch == "baseline":
                        continue
                    # Sanity check file meta start_ends
                    assert (
                        file_meta["meta"]["year_start"] == bin_range["file_year_start"]
                    )
                    assert file_meta["meta"]["year_end"] == bin_range["file_year_end"]
                    try:
                        tiff_fpath = self.generate_temporally_averaged_tiff(
                            file_meta,
                            epoch,
                            output_dir,
                            bin_range["bin_start"],
                            bin_range["bin_end"],
                        )
                    except Exception:
                        logging.exception("")
                    output_fpaths.append(tiff_fpath)
        return output_fpaths

    def occurrence_tiffs_to_exposure(
        self, occurrence_fpaths: List[str], output_dir: str
    ) -> List[str]:
        """
        Generate exposure tiffs from temporally averaged tiffs

        ::param occurrence_fpaths List[str] filepaths to occurrence tiffs
        """
        output_fpaths = []
        for occurrence_fpath in occurrence_fpaths:
            try:
                exposure_fpath = self.generate_popn_exposure_tiff(
                    occurrence_fpath, output_dir
                )
            except Exception:
                logging.exception("")
            logging.debug(
                "wrote exposure geotiff of size %s to: %s",
                os.path.getsize(exposure_fpath),
                exposure_fpath,
            )
            output_fpaths.append(exposure_fpath)
        return output_fpaths

    def run(
        self, download: bool = True, limit_files: int = None, log_meta: bool = False
    ) -> None:
        """
        Main runner - download files and process meta into CSV
        """
        logging.info("Running %s Downloader...", self.__class__.__name__)
        files_meta = self.fetch_nc_paths()
        if log_meta is True:
            logging.info(files_meta)
        if download is True:
            logging.info("downloading files...")
            files_meta = self.download_files(files_meta, limit_files=limit_files)
        # logging.info(("generating hazard csv")
        # h.append_hazard_csv(files_meta)
        print("Complete")


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    processor = HazardISIMPExtremeHeat(
        download_dir=args.download_dir,
        hazard_csv_fpath=args.hazard_csv_fpath,
        world_pop_2020_05deg_fpath=os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "worldpop_2020_1km_aggregated_resamp_0x5deg.tif",
        ),
    )
    # Generate the files metadata
    files_meta = processor.fetch_nc_paths()
    if args.log_meta == "True":
        logging.info(files_meta)
    # Download Files - checks local existance
    logging.info("downloading files...")
    result = processor.download_files(files_meta, limit_files=args.limit_files)
    if result is False:
        logging.warning("some files failed to download - see earlier messages")
    # generate averaged tiffs
    occurrence_tif_fpaths = processor.nc4_to_average_tiffs(
        files_meta, args.output_occurrence_data_directory
    )
    # Generate exposure tiffs
    exposure_tif_fpaths = processor.occurrence_tiffs_to_exposure(
        occurrence_tif_fpaths, args.output_exposure_data_directory
    )
    logging.info("Done.")
