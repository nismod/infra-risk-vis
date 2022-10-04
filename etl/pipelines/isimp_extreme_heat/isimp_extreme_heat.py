"""
Data Set Downloaders
"""

from ast import Assert
import sys
import os
import pprint
import csv
from html.parser import HTMLParser
from typing import List
from urllib.request import urlopen, urlretrieve
import logging
import argparse

logging.basicConfig(format="%(asctime)s %(levelname)s:%(message)s", level=logging.DEBUG)


parser = argparse.ArgumentParser(description="Downloaders")
parser.add_argument("dataset_type", type=str, help="Type of dataset to download")
parser.add_argument("hazard_csv_fpath", type=str, help="path to cumulative hazard file")
parser.add_argument("download_dir", type=str, help="Download directory for image files")
parser.add_argument(
    "--source_url", type=str, help="Base URL to use for Image file source"
)
parser.add_argument("--list_url", type=str, help="URL Listing files")
parser.add_argument(
    "--limit_files",
    type=int,
    default=None,
    help="Limit the number of files processed arbitrarily (for testing) - True / False",
)
parser.add_argument(
    "--download",
    default="True",
    type=str,
    help="Download the files (dry-run if False) - True / False",
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
        hazard_csv_fpath: str,
        download_dir: str,
        list_file: str = "https://github.com/nismod/open-gira/files/9488963/GRII_open_hazard_links.csv",
    ):
        self.list_file = list_file
        self.hazard_csv_fpath = hazard_csv_fpath
        self.download_dir = download_dir
        self.valid_hazard_types = ["extreme_heat"]
        raise Exception("incomplete")

    @staticmethod
    def nc_fname_from_url(url: str) -> str:
        """
        Parse filename for nc file from url
        """
        return url.split("/")[-1]

    def parse_fname(self) -> dict:
        """
        Parse the filename for hazard.csv keys - hazard, rp, rcp, epoch, gcm
        """

    def fetch_nc_paths(self) -> List[dict]:
        """
        Fetch the NC CSV and parse out the ISIMP Extreme Heat paths
        """
        # fetch list to download dir then remove
        local_list_file = os.path.join(self.download_dir, "tmp_list")
        try:
            fsize = download_file(self.list_file, local_list_file)
            if not fsize:
                raise FileNotFoundError("failed to download list file")
            with open(local_list_file, "r") as csvfile:
                reader = csv.DictReader(csvfile)
                files_meta = list(reader)
                # # Remove extraneous types
                files_meta = [
                    item
                    for item in files_meta
                    if item["\ufeffhazard"] in self.valid_hazard_types
                ]
                # Append the filenames
                for item in files_meta:
                    item["fname"] = self.nc_fname_from_url(item["url"])
                return files_meta
        except:
            raise
        finally:
            os.remove(local_list_file)

    def download_files(self, files_meta: List[dict], limit_files=None) -> bool:
        """
        Download the files in the given meta

        ::returns bool All completed
        """
        fsizes = []
        for idx, file in enumerate(files_meta):
            logging.debug("Downloading %s", file["url"])
            try:
                fpath = os.path.join(self.download_dir, file["fname"])
                # Skip if file already exists
                if os.path.exists(fpath):
                    logging.info(
                        "file already exists locally - skipping download: %s",
                        file["fname"],
                    )
                if limit_files and idx >= limit_files:
                    logging.info("download aborting due to limit_files=%s", limit_files)
                    break

                fsize = download_file(file["url"], fpath)
                if not fsize:
                    raise Exception(
                        "download resulted in invalid file: {}".format(file["url"])
                    )
                fsizes.append(fsize)
            except Exception as err:
                logging.exception("")
        if limit_files:
            return len(fsizes) == limit_files
        return len(fsizes) == len(files_meta)

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
            result = self.download_files(files_meta, limit_files=limit_files)
            if result is False:
                logging.warning("some files failed to download - see earlier messages")
        # logging.info(("generating hazard csv")
        # h.append_hazard_csv(files_meta)
        print("Complete")


if __name__ == "__main__":
    args = parser.parse_args()
    if args.dataset_type == "aqueduct":
        h = HazardAqueduct(args.hazard_csv_fpath, args.download_dir)
        # fnames = h.fetch_tiff_fnames()
        # files_meta = h.parse_fnames(fnames)
        # pp = pprint.PrettyPrinter(indent=4)
        # pp.pprint(files_meta)
        # # h.append_hazard_csv(files_meta)
        h.run(
            download=args.download == "True",
            limit_files=args.limit_files,
            log_meta=args.log_meta == "True",
        )
    elif args.dataset_type == "isimp_extreme_heat":
        h = HazardISIMPExtremeHeat(args.hazard_csv_fpath, args.download_dir)
        h.run(
            download=args.download == "True",
            limit_files=args.limit_files,
            log_meta=args.log_meta == "True",
        )
    else:
        logging.error("Unknown Downloader type:", args.dataset_type)
