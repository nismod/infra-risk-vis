"""
Aqueduct dataset parser, takes remote listing and creates CSV file with rows of
rasters and their metadata.
"""

import sys
import os
import pprint
import csv
from html.parser import HTMLParser
from typing import List, Tuple
from urllib.request import urlopen, urlretrieve
import argparse

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
)

from pipelines.helpers import get_logger

LOG = get_logger(__name__)


parser = argparse.ArgumentParser(description="Aqueduct Downloader")
parser.add_argument(
    "--source_url", type=str, help="Base URL to use for Image file source"
)
parser.add_argument("--list_url", type=str, help="URL Listing files")
parser.add_argument(
    "--log_meta",
    type=bool,
    default=False,
    help="Logout metadata associated with files",
)
parser.add_argument(
    "--csv_path",
    type=str,
    default=None,
    help="File-path at-which to dump file meta to CSV",
)


class ListHTMLParser(HTMLParser):
    files = []
    accepted_extensions = ["tif"]

    def handle_data(self, data):
        try:
            if data.split(".")[1] in self.accepted_extensions:
                self.files.append(data)
        except Exception as err:
            pass


class HazardAqueduct:
    """
    Pre-requisites:
        - Hosting URL for Aqueduct data (Costal and Riverine)
    """

    def __init__(
        self,
        source_url: str = "http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2",
        list_url: str = "http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/index.html",
        hazard_csv_fpath: str = None,
        skip_fname_patterns: List = [
            "perc_05",
            "0_perc_05",
            "perc_50",
            "0_perc_50",
            "nosub",
        ],
    ):
        self.source_url = source_url
        self.list_url = list_url
        self.hazard_csv_fpath = hazard_csv_fpath
        self.hazard_csv_fieldnames = [
            "hazard",
            "path",
            "rp",
            "rcp",
            "epoch",
            "gcm",
            "key",
            "url",
        ]
        # If these patterns are in the filenames they will be removed from processing
        self.skip_fname_patterns = skip_fname_patterns
        # Top level ETL direcctory
        self.etl_dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    def coastal_fname_map(self):
        """
        Filename map for coastal data to CSV variables
        See: https://onedrive.live.com/edit.aspx?resid=55B6ADF51ADBA854!193&cid=55b6adf51adba854&CT=1664373823410&OR=ItemsView
        """
        # idx: {input : output}
        return {
            "hazard": "coastal",
            "fname_map": {
                1: {
                    "input": "climatescenario",
                    "output": "rcp",
                    "map": {"historical": "baseline", "rcp4p5": "4x5", "rcp8p5": "8x5"},
                },
                3: {
                    "input": "year",
                    "output": "epoch",
                    "map": {
                        "hist": "present",
                        "1980": "present",
                        "2030": "2030",
                        "2050": "2050",
                        "2080": "2080",
                    },
                },
                4: {
                    "input": "returnperiod",
                    "output": "rp",
                    "map": {
                        "rp0001": "1",
                        "rp00001": "1",
                        "rp0002": "2",
                        "rp00002": "2",
                        "rp0005": "5",
                        "rp00005": "5",
                        "rp0010": "10",
                        "rp00010": "10",
                        "rp0025": "25",
                        "rp00025": "25",
                        "rp0050": "50",
                        "rp00050": "50",
                        "rp0100": "100",
                        "rp00100": "100",
                        "rp0250": "250",
                        "rp00250": "250",
                        "rp0500": "500",
                        "rp00500": "500",
                        "rp1000": "1000",
                        "rp01000": "1000",
                    },
                },
            },
        }

    def fluvial_fname_map(self):
        """
        Filename map for fluvial data to CSV variables
        See: http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/_Aqueduct_Floods_Data_Dictionary.xlsx
        """
        # idx: {input : output}
        return {
            "hazard": "fluvial",
            "fname_map": {
                1: {
                    "input": "climatescenario",
                    "output": "rcp",
                    "map": {"historical": "baseline", "rcp4p5": "4x5", "rcp8p5": "8x5"},
                },
                2: {
                    "input": "model",
                    "output": "gcm",
                    "map": {
                        "000000000WATCH": "WATCH",
                        "00000NorESM1-M": "NorESM1-M",
                        "0000GFDL_ESM2M": "GFDL_ESM2M",
                        "0000GFDL-ESM2M": "GFDL-ESM2M",
                        "0000HadGEM2-ES": "HadGEM2-ES",
                        "00IPSL-CM5A-LR": "IPSL-CM5A-LR",
                        "MIROC-ESM-CHEM": "MIROC-ESM-CHEM",
                    },
                },
                3: {
                    "input": "year",
                    "output": "epoch",
                    "map": {
                        "hist": "present",
                        "1980": "present",
                        "2030": "2030",
                        "2050": "2050",
                        "2080": "2080",
                    },
                },
                4: {
                    "input": "returnperiod",
                    "output": "rp",
                    "map": {
                        "rp0001": "1",
                        "rp00001": "1",
                        "rp0002": "2",
                        "rp00002": "2",
                        "rp0005": "5",
                        "rp00005": "5",
                        "rp0010": "10",
                        "rp00010": "10",
                        "rp0025": "25",
                        "rp00025": "25",
                        "rp0050": "50",
                        "rp00050": "50",
                        "rp0100": "100",
                        "rp00100": "100",
                        "rp0250": "250",
                        "rp00250": "250",
                        "rp0500": "500",
                        "rp00500": "500",
                        "rp1000": "1000",
                        "rp01000": "1000",
                    },
                },
            },
        }

    def build_tiff_link(self, fname: str) -> str:
        """
        Build a link to given tiff
        """
        return "{}/{}".format(self.source_url, fname)

    def fetch_tiff_fnames(self) -> List[str]:
        """
        Extract all tiff links from the source URL
        """
        parser = ListHTMLParser()
        with urlopen(self.list_url) as f:
            content = f.read()
            # Extract the fnames from source
            parser.feed(content.decode())
        # Generate full links
        return parser.files

    def map_selector(self, fname: str):
        """
        Select which mapper to use forthe given fname
        """
        _type = fname.split("_")[0]
        try:
            return {
                "inuncoast": self.coastal_fname_map(),
                "inunriver": self.fluvial_fname_map(),
            }[_type]
        except KeyError:
            raise Exception("unknown map type: {}".format(_type))

    def parse_fname(self, fname: str) -> dict:
        """
        Parse a fname through the relevent mapper

        ::return meta dict
            {   'filename': 'inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif',
                'meta': {   'epoch': '2080',
                            'gcm': 'MIROC-ESM-CHEM',
                            'hazard': 'fluvial',
                            'key': 'inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000',
                            'path': 'input/hazard-aqueduct/global/inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif',
                            'rcp': '8.5'},
                'url': 'http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif'
            }
        """
        data_mapper = self.map_selector(fname)
        output = {
            "filename": fname,
            "meta": {
                "url": self.build_tiff_link(fname),
            }
        }
        try:
            for idx, item in enumerate(os.path.splitext(fname)[0].split("_")):
                if idx == 2 and data_mapper["hazard"] == "coastal":
                    # No GCM for Coastal
                    output["meta"]["gcm"] = "None"
                elif idx in data_mapper["fname_map"].keys():
                    output["meta"][
                        data_mapper["fname_map"][idx]["output"]
                    ] = data_mapper["fname_map"][idx]["map"][item]
                    # Add hazard type
                    output["meta"]["hazard"] = data_mapper["hazard"]
                else:
                    continue
            if not output:
                raise Exception("meta failed to parse any indices for file")
        except Exception as err:
            print("failed to parse file", fname, idx, item, err)
        return output

    def parse_fnames(self, fnames: List[str]) -> List[dict]:
        """
        Parse the filenames into mapped elements
        """
        output = []
        for fname in fnames:
            try:
                # Skip if pattern matches skip cases
                if any(checkstr in fname for checkstr in self.skip_fname_patterns):
                    print(f"Skipping {fname} due to substring in skip-list")
                    continue
                # Generate the meta from filename
                file_meta = self.parse_fname(fname)
                file_meta["meta"]["key"] = os.path.splitext(fname)[0]
            except Exception as err:
                print("failed to parse file:", fname, err)
            output.append(file_meta)
        return output

    def file_metadata(self) -> List[dict]:
        """
        Generate metadata for source files

        ::returns List[dict] metadata for files:
            {
                'filename': 'inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif',
                'url': 'http://wri-projects.s3.amazonaws.com/AqueductFloodTool/download/v2/inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000.tif',
                'meta': {
                    'rcp': '8x5',
                    'hazard': 'fluvial',
                    'gcm': 'MIROC-ESM-CHEM',
                    'epoch': '2080',
                    'rp': '1000',
                    'key': 'inunriver_rcp8p5_MIROC-ESM-CHEM_2080_rp01000'
                }
            }
        """
        fnames = self.fetch_tiff_fnames()
        LOG.info("parsing file meta")
        files_meta = self.parse_fnames(fnames)
        return files_meta

    def hazard_csv_exists(self) -> bool:
        """
        Check the provided hazard csv exists
        """
        return os.path.exists(self.hazard_csv_fpath)

    def hazard_csv_valid(self) -> bool:
        """
        Ensure the provided hazard CSV is valid
        """
        with open(self.hazard_csv_fpath, "r") as csvfile:
            d_reader = csv.DictReader(csvfile)
            headers = d_reader.fieldnames
            assert set(headers) == set(
                self.hazard_csv_fieldnames
            ), "existing hazard csv has header mismatch: {} | {}".format(
                headers, self.hazard_csv_fieldnames
            )
            return True

    def count_hazard_csv_rows(self) -> int:
        """
        Count number of rows in the Hazard CSV
        """
        return sum(1 for _ in open(self.hazard_csv_fpath))

    def write_hazard_csv(self, files_meta: List[dict]) -> None:
        """
        Write the meta from parsed raster list to the hazard CSV table
        """
        # overwrite/create file if req
        with open(self.hazard_csv_fpath, "w") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self.hazard_csv_fieldnames)
            writer.writeheader()

            # Dump meta
            for file in files_meta:
                row = file["meta"]
                row["path"] = file["filename"]
                writer.writerow(row)

        # Check CSV is valid (headers etc)
        _ = self.hazard_csv_valid()

        csv_count = self.count_hazard_csv_rows()
        print(
            "CSV length:",
            csv_count,
            "Number files processed:",
            len(files_meta),
        )

    def run(
        self,
        download_path: str = None,
        log_meta: bool = False,
        csv_path: str = None,
        **kwargs
    ) -> None:
        """
        Main runner - download files and process meta into CSV
        """
        print("Running HazardAqueduct Downloader...")
        print("fetch file meta")
        fnames = self.fetch_tiff_fnames()
        print("parse file meta")
        files_meta = self.parse_fnames(fnames)
        if log_meta:
            pp = pprint.PrettyPrinter(indent=4)
            pp.pprint(files_meta)
        if download_path:
            print("downloading files...")
            files_meta = self.download_files(files_meta)
        if csv_path:
            print("generating hazard csv")
            self.write_hazard_csv(files_meta)
        print("Complete")


if __name__ == "__main__":
    parsed_args = vars(parser.parse_args())
    processor = HazardAqueduct(
        hazard_csv_fpath=parsed_args["csv_path"],
        source_url=parsed_args["source_url"],
        list_url=parsed_args["list_url"]     
    )
    processor.run(**parsed_args)