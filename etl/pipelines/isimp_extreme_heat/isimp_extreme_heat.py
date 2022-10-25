"""
Data Set Downloaders
"""

import os
import sys
import csv
from typing import List
import logging
import argparse
from dataclasses import dataclass, asdict

import numpy as np
from netCDF4 import Dataset
from osgeo import gdal, osr, gdalconst

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from pipelines.helpers import (
    count_hazard_csv_rows,
    download_file,
    nc_fname_from_url,
    fname_from_fpath,
    file_key_from_fname,
    hazard_csv_exists,
    hazard_csv_valid,
)

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
    "--log_meta",
    type=str,
    default="False",
    help="Logout metadata associated with files",
)


@dataclass
class ISIMPExtremeHeatMeta:
    """Metadata associated with an input file"""

    hazard: str
    key: str
    fname: str
    maintainer: str
    model: str
    climate_forcing: str
    bias_adj: str
    climate_scenario: str
    socio_econ_scenario: str
    sensitivity_scenario: str
    variable: str
    region: str
    time_step: str
    year_start: int
    year_end: int

    def external_rcp(self):
        """
        'rcp60' -> 6x0
        """
        if self.climate_scenario == "historical":
            return "baseline"
        if "rcp" not in self.climate_scenario:
            raise Exception(f"rcp field appears to be wrong for {self.fname}")
        base = self.climate_scenario.replace("rcp", "")
        return f"{base[0]}x{base[1]}"


@dataclass
class ISIMPExtremeHeatOccurrence:
    """Metadata associated with a processed occurrence file"""

    metric = "occurrence"
    source_fname: str
    fname: str
    fpath: str
    key: str
    epoch: int
    gcm: str
    rcp: str


@dataclass
class ISIMPExtremeHeatExposure:
    """Metadata associated with a processed exposure file"""

    metric = "exposure"
    source_fname: str
    fname: str
    fpath: str
    key: str
    epoch: int
    gcm: str
    rcp: str


class ISIMPExtremeHeatFile:
    def __init__(self, url: str):
        self.source = "ISIMP"
        self.url = url
        self.fpath = None
        self.fsize = None
        self.fname = nc_fname_from_url(self.url)
        self.meta = self._parse_fname(self.fname)
        # Outputs from the occurrence and exposure processing related to this file
        self.occurrence_outputs: List[ISIMPExtremeHeatOccurrence] = []
        self.exposure_outputs: List[ISIMPExtremeHeatExposure] = []

    def _parse_fname(self, fname: str) -> ISIMPExtremeHeatMeta:
        """
        Parse the fname for hazard.csv keys

        ::param fname st e.g. lange2020_hwmid-humidex_gfdl-esm2m_ewembi_historical_nosoc_co2_leh_global_annual_1861_2005.nc4
            {MAINTAINER}_{MODEL}_{CLIMATE_FORCING}_{BIAS_ADJ}_{CLIMATE_SCENARIO(RCP)}_{SCO_ECO_SCENARIO}_{SENS_SCENARIO}_{VARIABLE}_{REGION}_{TIME_STEP}_{YEAR_START}_{YEAR_END}.nc4

        ::returns meta ISIMPExtremeHeatMeta
        """
        _key = file_key_from_fname(fname)
        parts = _key.split("_")
        return ISIMPExtremeHeatMeta(
            hazard="extreme_heat",
            fname=fname,
            key=_key,
            maintainer=parts[0],
            model=parts[1],
            climate_forcing=parts[2],
            bias_adj=parts[3],
            climate_scenario=parts[4],
            socio_econ_scenario=parts[5],
            sensitivity_scenario=parts[6],
            variable=parts[7],
            region=parts[8],
            time_step=parts[9],
            year_start=int(parts[10]),
            year_end=int(parts[11]),
        )


class HazardISIMPExtremeHeat:
    """
    ISIMP Extreme Heat
    """

    def __init__(
        self,
        download_dir: str = None,
        hazard_csv_fpath: str = None,
        world_pop_fpath: str = "GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0_wgs84.tif",
    ):
        self.hazard_csv_fpath = hazard_csv_fpath
        self.download_dir = download_dir
        self.world_pop_fpath = world_pop_fpath
        self.hazard_name = "extreme_heat"
        self.tmp_dir = "/tmp"  # For list file processing
        self.hazard_csv_fieldnames = [
            "hazard",
            "metric",  # occurrence or exposure
            "path",
            "rcp",
            "epoch",
            "gcm",
            "key",
        ]
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
        self.rcp_to_scenario_map = [
            "baseline",
            "2x6",
            "6x0",
        ]  # these map to indices of climate_scenarios_to_process
        self.time_spans_to_process = [[1861, 2005], [2006, 2099]]
        if not os.path.exists(self.world_pop_fpath):
            logging.warning(
                "world pop data does not exist at given path - exposure calcs will fail: %s",
                self.world_pop_fpath,
            )

    def fetch_nc_paths(self, remote_list_file: str) -> List[ISIMPExtremeHeatFile]:
        """
        Fetch the NC CSV and parse out the ISIMP Extreme Heat paths and file meta

        ::param remote_list_file str The remote list file to fetch

        ::returns file meta List[dict] e.g.
            [
                {'\ufeffhazard': 'extreme_heat',
                'source': 'ISIMIP',
                'url': 'https://files.isimip.org/ISIMIP2b/DerivedOutputData/Lange2020/HWMId-humidex/gfdl-esm2m/future/lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099.nc4',
                'fname': 'lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099.nc4',
                'path': '',
                'meta': {
                    'hazard': 'extreme_heat',
                    'key': 'lange2020_hwmid-humidex_gfdl-esm2m_ewembi_rcp26_nosoc_co2_leh_global_annual_2006_2099',
                    'maintainer': 'lange2020',
                    'model': 'hwmid-humidex',
                    'climate_forcing': 'gfdl-esm2m',
                    'bias_adj': 'ewembi',
                    'climate_scenario': 'rcp26',
                    'socio_econ_scenario': 'nosoc',
                    'sensitivity_scenario': 'co2',
                    'variable': 'leh',
                    'region': 'global',
                    'time_step': 'annual',
                    'year_start': 2006,
                    'year_end': 2099}
                }
            ]
        """
        input_files = []
        # fetch list to download dir then remove
        local_list_file = os.path.join(self.tmp_dir, "extreme_heat_list.csv")
        try:
            fsize = download_file(remote_list_file, local_list_file)
            if not fsize:
                raise FileNotFoundError("failed to download list file")
            with open(local_list_file, "r") as csvfile:
                reader = csv.DictReader(csvfile)
                files = list(reader)
                # Remove extraneous types, Append the filenames and parse filename to meta
                for item in files:
                    if item["\ufeffhazard"] == self.hazard_name:
                        input_file = ISIMPExtremeHeatFile(item["url"])
                        # Remove climate scenarios and years we do not require
                        if (
                            input_file.meta.climate_scenario
                            in self.climate_scenarios_to_process
                        ):
                            if [
                                input_file.meta.year_start,
                                input_file.meta.year_end,
                            ] in self.time_spans_to_process:
                                input_files.append(input_file)
                return input_files
        except:
            raise
        finally:
            os.remove(local_list_file)

    def download_files(
        self, input_files: List[ISIMPExtremeHeatFile]
    ) -> List[ISIMPExtremeHeatFile]:
        """
        Download the files in the given meta

        ::returns input_files - updated with local file paths
        """
        for idx, input_file in enumerate(input_files):
            logging.debug("Downloading %s", input_file.url)
            try:
                fpath = os.path.join(self.download_dir, input_file.fname)
                # Skip if file already exists
                if os.path.exists(fpath):
                    logging.info(
                        "file already exists locally - skipping download: %s",
                        input_file.fname,
                    )
                    fsize = os.path.getsize(fpath)
                else:
                    fsize = download_file(input_file.url, fpath)
                    if not fsize:
                        raise Exception(
                            "download resulted in invalid file: {}".format(
                                input_file.url
                            )
                        )
                input_file.fsize = fsize
                input_file.fpath = fpath
            except Exception:
                logging.exception("")
        return input_files

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
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        ds = drv.Create(output_fpath, width, height, 1, gdal.GDT_Float32)
        ds.SetProjection(srs.ExportToWkt())
        ds.SetGeoTransform((ul_x, pixel_size, 0, ul_y, 0, -pixel_size))
        ds.GetRasterBand(1).WriteArray(data)

    def occurrence_tiff_fname_from_meta(
        self, input_file: ISIMPExtremeHeatFile, epoch: str
    ) -> str:
        """
        Generate a name for a processed tiff using the files metadata
        """
        basename = file_key_from_fname(input_file.fname)
        return f"{basename}_{epoch}_occurrence.tif"

    def exposure_tiff_fname_from_occurrence(self, occurrence_tif_fname: str) -> str:
        """
        Generate a name for a exposure tif from occurrence
        """
        basename = file_key_from_fname(occurrence_tif_fname)
        fname = basename.replace("occurrence", "exposure")
        return f"{fname}.tif"

    def match_raster_resoluton_and_extent(
        self, master_fpath: str, slave_fpath: str, output_fpath: str, res=None
    ) -> bool:
        """
        Generate a new raster that is a copy of slave, but with the resolution and size of the master
        """
        slave_ds = gdal.Open(slave_fpath)
        if slave_ds is None:
            raise IOError
        slave_proj = slave_ds.GetProjection()
        data_type = slave_ds.GetRasterBand(1).DataType
        n_bands = slave_ds.RasterCount

        master_ds = gdal.Open(master_fpath)
        if master_ds is None:
            raise IOError
        master_proj = master_ds.GetProjection()
        master_geotrans = master_ds.GetGeoTransform()
        w = master_ds.RasterXSize
        h = master_ds.RasterYSize
        if res is not None:
            master_geotrans[1] = float(res)
            master_geotrans[-1] = -float(res)

        dst_ds = gdal.GetDriverByName("GTiff").Create(
            output_fpath, w, h, n_bands, data_type, ["COMPRESS=LZW"]
        )
        dst_ds.SetGeoTransform(master_geotrans)
        dst_ds.SetProjection(master_proj)

        gdal.ReprojectImage(
            slave_ds, dst_ds, slave_proj, master_proj, gdal.GRA_NearestNeighbour
        )
        dst_ds = None  # Flush to disk
        return os.path.exists(output_fpath)

    def generate_popn_exposure_tiff(
        self, input_occurrence_file: ISIMPExtremeHeatOccurrence, output_dir: str
    ) -> ISIMPExtremeHeatExposure:
        """
        Generate a TIFF file of world population x heatwave occurrence
        """
        logging.info(
            "generating popn exposure tiff - input file: %s",
            input_occurrence_file.fpath,
        )
        # Generate the output filename
        output_fname = self.exposure_tiff_fname_from_occurrence(
            os.path.basename(input_occurrence_file.fpath)
        )
        output_fpath = os.path.join(output_dir, output_fname)
        # Skip if it exists
        if os.path.exists(output_fpath):
            logging.warning(
                "pop exposure output fpath appears to exist - skipping reprocessing: %s",
                output_fpath,
            )
            return ISIMPExtremeHeatExposure(
                input_file.fname,
                output_fname,
                output_fpath,
                file_key_from_fname(output_fname),
                input_occurrence_file.epoch,
                input_occurrence_file.gcm,
                input_occurrence_file.rcp,
            )

        logging.debug(
            "generating popn exposure tiff - opening worldpop: %s",
            self.world_pop_fpath,
        )

        # Warp Occurrence to Worldpop proj, cellsize and extent
        logging.debug(
            "generating popn exposure tiff - warping occurrence to world pop size and proj"
        )
        tmp_exposure_fpath = os.path.join(output_dir, "tmp_reproj.tif")
        result = self.match_raster_resoluton_and_extent(
            self.world_pop_fpath, input_occurrence_file.fpath, tmp_exposure_fpath
        )
        logging.debug(
            "generating popn exposure tiff - matched res and scale of world pop raster: %s",
            result,
        )

        logging.debug("generating popn exposure tiff - running exposure calc...")
        # Run the calculation for exposure
        calc_cmd = f'gdal_calc.py -A {self.world_pop_fpath} -B {tmp_exposure_fpath} --co="COMPRESS=LZW" --outfile={output_fpath} --calc="A*B"'
        os.system(calc_cmd)

        logging.debug("generating popn exposure tiff - removing tmp files")
        os.remove(tmp_exposure_fpath)

        return ISIMPExtremeHeatExposure(
            input_file.fname,
            output_fname,
            output_fpath,
            file_key_from_fname(output_fname),
            input_occurrence_file.epoch,
            input_occurrence_file.gcm,
            input_occurrence_file.rcp,
        )

    def generate_temporally_averaged_tiff(
        self,
        input_file: ISIMPExtremeHeatFile,
        output_epoch: str,
        output_dir: str,
        averaging_year_start: int,
        averaging_year_end: int,
    ) -> ISIMPExtremeHeatOccurrence:
        """
        Generate a TIFF file using temporal averaging for the given time window
        """
        logging.info(
            "generating temporally averaged tiff - input file: %s epoch: %s: averaging_year_start: %s, averaging_year_end: %s",
            input_file.fname,
            output_epoch,
            averaging_year_start,
            averaging_year_end,
        )
        output_fname = self.occurrence_tiff_fname_from_meta(input_file, output_epoch)
        output_fpath = os.path.join(output_dir, output_fname)
        # Load nc4
        dataset = Dataset(input_file.fpath, "r", format="NETCDF4")
        # Load data
        leh_np = dataset.variables["leh"][:]
        # Calculate temporal slice idx
        start_idx = averaging_year_start - input_file.meta.year_start
        end_idx = averaging_year_end - input_file.meta.year_start
        logging.debug(
            "slicing file %s at idx %s (yr %s) to %s (yr %s)",
            input_file.fname,
            start_idx,
            input_file.meta.year_start + start_idx,
            end_idx,
            input_file.meta.year_start + end_idx,
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
        logging.debug("writing occurrence geotiff to: %s", output_fpath)
        self.np_to_geotiff(leh_np_mean, output_fpath)
        logging.debug(
            "wrote occurrence geotiff of size %s to: %s",
            os.path.getsize(output_fpath),
            output_fpath,
        )
        return ISIMPExtremeHeatOccurrence(
            input_file.fname,
            output_fname,
            output_fpath,
            file_key_from_fname(output_fname),
            output_epoch,
            input_file.meta.climate_forcing,
            input_file.meta.external_rcp(),
        )

    def nc4_to_average_tiffs(
        self, input_files: List[ISIMPExtremeHeatFile], output_dir: str
    ) -> List[ISIMPExtremeHeatFile]:
        """
        Temporal averaging for netCDF4 files contained in files meta
            Ref: Dr Mark bernhofen

        Note we only process the given climate senario files

        ::param input_files List[ISIMPExtremeHeatFile] - metadata generated from list_url files
        ::param output_dir str - Processing output directory

        ::return input_files List[ISIMPExtremeHeatFile] with occurrence datasets populated
        """
        for input_file in input_files:
            # In case some extraneous files are present:
            if (
                input_file.meta.climate_scenario
                not in self.climate_scenarios_to_process
            ):
                continue
            if input_file.meta.climate_scenario == "historical":
                # Sanity check file meta start_ends
                assert (
                    input_file.meta.year_start
                    == self.epoch_bins["baseline"]["file_year_start"]
                )
                assert (
                    input_file.meta.year_end
                    == self.epoch_bins["baseline"]["file_year_end"]
                )
                try:
                    # Special case for historical
                    occurrence_meta = self.generate_temporally_averaged_tiff(
                        input_file,
                        "baseline",
                        output_dir,
                        self.epoch_bins["baseline"]["bin_start"],
                        self.epoch_bins["baseline"]["bin_end"],
                    )
                    input_file.occurrence_outputs.append(occurrence_meta)
                except Exception:
                    logging.exception("")
            else:
                for epoch, bin_range in self.epoch_bins.items():
                    if epoch == "baseline":
                        continue
                    # Sanity check file meta start_ends
                    assert input_file.meta.year_start == bin_range["file_year_start"]
                    assert input_file.meta.year_end == bin_range["file_year_end"]
                    try:
                        occurrence_meta = self.generate_temporally_averaged_tiff(
                            input_file,
                            epoch,
                            output_dir,
                            bin_range["bin_start"],
                            bin_range["bin_end"],
                        )
                        input_file.occurrence_outputs.append(occurrence_meta)
                    except Exception:
                        logging.exception("")
        return input_files

    def occurrence_tiffs_to_exposure(
        self, input_files: List[ISIMPExtremeHeatFile], output_dir: str
    ) -> List[ISIMPExtremeHeatFile]:
        """
        Generate exposure tiffs from temporally averaged tiffs

        ::param input_files List[ISIMPExtremeHeatFile] ISIMPExtremeHeatFile objects containing meta for occurrence datasets
        """
        for input_file in input_files:
            for occurrence_file in input_file.occurrence_outputs:
                try:
                    exposure = self.generate_popn_exposure_tiff(
                        occurrence_file, output_dir
                    )
                    input_file.exposure_outputs.append(exposure)
                    logging.debug(
                        "wrote exposure geotiff of size %s to: %s",
                        os.path.getsize(exposure.fpath),
                        exposure.fpath,
                    )
                except Exception:
                    logging.exception("")
        return input_files

    def append_hazard_csv(self, input_files: List[ISIMPExtremeHeatFile]) -> None:
        """
        Append the meta from parsed files to the hazard csv
        """
        # Generate file if req
        if not hazard_csv_exists(self.hazard_csv_fpath):
            with open(self.hazard_csv_fpath, "w") as csvfile:
                print(" No hazard csv found writing new at: ", self.hazard_csv_fpath)
                # Generate a new file
                writer = csv.DictWriter(csvfile, fieldnames=self.hazard_csv_fieldnames)
                writer.writeheader()
        # Check CSV is valid (headers etc)
        _ = hazard_csv_valid(self.hazard_csv_fpath, self.hazard_csv_fieldnames)
        # Dump meta
        count_before = count_hazard_csv_rows(self.hazard_csv_fpath)
        with open(self.hazard_csv_fpath, "a") as csvfile:
            # Setup the writer
            writer = csv.DictWriter(csvfile, fieldnames=self.hazard_csv_fieldnames)
            # Dump meta
            count_outputs = 0
            for input_file in input_files:
                row = {}
                row["hazard"] = input_file.meta.hazard
                for occurrence_file in input_file.occurrence_outputs:
                    # Occurrence
                    row["rcp"] = occurrence_file.rcp
                    row["gcm"] = occurrence_file.gcm
                    row["epoch"] = occurrence_file.epoch
                    row["key"] = occurrence_file.key
                    row["metric"] = occurrence_file.metric
                    row["path"] = occurrence_file.fpath
                    writer.writerow(row)
                    count_outputs += 1
                # Exposure
                for exposure_file in input_file.exposure_outputs:
                    row["rcp"] = exposure_file.rcp
                    row["gcm"] = exposure_file.gcm
                    row["epoch"] = exposure_file.epoch
                    row["key"] = exposure_file.key
                    row["metric"] = exposure_file.metric
                    row["path"] = exposure_file.fpath
                    writer.writerow(row)
                    count_outputs += 1
        count_after = count_hazard_csv_rows(self.hazard_csv_fpath)
        print(
            "Count before:",
            count_before - 1,
            "Count after:",
            count_after,
            "number input files processed:",
            len(input_files),
        )


if __name__ == "__main__":
    args = parser.parse_args()
    print(args)
    processor = HazardISIMPExtremeHeat(
        download_dir=args.download_dir,
        hazard_csv_fpath=args.hazard_csv_fpath,
        world_pop_fpath=os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "common",
            "data",
            "GHS_POP_E2020_GLOBE_R2022A_54009_1000_V1_0_wgs84.tif",
        ),
    )
    # Generate the files metadata
    input_files = processor.fetch_nc_paths(
        "https://github.com/nismod/open-gira/files/9488963/GRII_open_hazard_links.csv"
    )
    if args.log_meta == "True":
        for input_file in input_files:
            logging.info(asdict(input_file.meta))
    # Download Files - checks local existance
    logging.info("downloading files...")
    input_files = processor.download_files(input_files)
    # generate averaged tiffs
    input_files = processor.nc4_to_average_tiffs(
        input_files, args.output_occurrence_data_directory
    )
    for input_file in input_files:
        print(
            input_file.fname,
            len(input_file.occurrence_outputs),
            len(input_file.exposure_outputs),
        )
    # Generate exposure tiffs
    input_files = processor.occurrence_tiffs_to_exposure(
        input_files, args.output_exposure_data_directory
    )
    # Generate the CSV
    processor.append_hazard_csv(input_files)
    logging.info("Done.")
