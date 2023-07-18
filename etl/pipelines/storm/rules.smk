rule download_raw_data:
    """
    Download files from remote location.
    """
    output:
        # requesting any file will trigger a download of the whole archive
        requested_raster = "raster/raw/storm/{KEY}.tif",
    shell:
        """
        zenodo_get --record=7438145 --output-dir=$(dirname {output.requested_raster})
        """