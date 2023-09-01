"""
Download resources from the Copernicus CDS API.
"""

import argparse
import os

import cdsapi


def download_from_CDS(dataset_name: str, variable: str, file_format: str, version: str, year: str, output_path: str) -> None:
    """
    Download a resource from the Copernicus CDS API, given appropriate credentials.

    Requires COPERNICUS_CDS_URL and COPERNICUS_CDS_API_KEY to be in the environment.
    For more details see: https://cds.climate.copernicus.eu/api-how-to

    Args:
        dataset_name: Name of dataset to download
        variable: Name of variable to request
        file_format: Desired file format e.g. zip
        version: Version of dataset
        year: Year of dataset applicability
        output_path: Where to save the downloaded file
    """

    client = cdsapi.Client(
        url=os.environ.get("COPERNICUS_CDS_URL"),
        key=os.environ.get("COPERNICUS_CDS_API_KEY")
    )

    # N.B. Files are covered by licences which need to be manually accepted, e.g.
    # https://cds.climate.copernicus.eu/cdsapp/#!/terms/satellite-land-cover
    # https://cds.climate.copernicus.eu/cdsapp/#!/terms/vito-proba-v
    #
    # Ideally we could programmatically accept the necessary licence conditions
    # the below code is an attempt at that, but fails with an HTTP 403, not
    # logged in when trying to simulate a user acceptance 
    #
    #   API_URL = os.environ.get("COPERNICUS_CDS_URL")
    #   payloads = [
    #       [{"terms_id":"vito-proba-v","revision":1}],
    #       [{"terms_id":"satellite-land-cover","revision":1}],
    #   ]
    #   for payload in payloads:
    #       client._api(
    #           url=f"{API_URL.rstrip('/')}.ui/user/me/terms-and-conditions",
    #           request=payload,
    #           method="post"
    #       )
    #
    # See https://github.com/ecmwf/cdsapi/blob/master/cdsapi/api.py

    client.retrieve(
        dataset_name,
        {
            'variable': variable,
            'format': file_format,
            'version': version,
            'year': year,
        },
        output_path
    )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--dataset_name", help="Name of dataset to download")
    parser.add_argument("-V", "--variable", default="all", help="Variable of dataset to download")
    parser.add_argument("-f", "--file_format", default="zip", help="File format to download dataset in")
    parser.add_argument("-v", "--version", help="Version of dataset to download")
    parser.add_argument("-y", "--year", help="Applicable year of dataset")
    parser.add_argument("-o", "--output_path", help="File path to write dataset to")

    kwargs = vars(parser.parse_args())
    download_from_CDS(**kwargs)