rule download_raw_data:
    """
    Download files from remote location.
    """
    output:
        # requesting any file will trigger a download of the whole archive
        requested_raster = "raster/raw/storm/{KEY}.tif",
    shell:
        """
        OUTPUT_DIR=$(dirname {output.requested_raster})

        pushd $OUTPUT_DIR
            zenodo_get -w links.txt --record=7438145
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """
