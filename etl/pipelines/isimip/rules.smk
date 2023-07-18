rule download:
    """
    Fetch our preprocessed data from zenodo.

    N.B. Requesting any file will trigger a download of the whole archive
    """
    output:
        "raster/raw/isimip/{KEY}.tif",
    shell:
        """
        OUTPUT_DIR=$(dirname {output})
        ARCHIVE=$OUTPUT_DIR/lange2020_expected_occurrence.zip

        zenodo_get --record=8147088 --output-dir=$OUTPUT_DIR

        unzip -d $OUTPUT_DIR

        mv $OUTPUT_DIR/data/* $OUTPUT_DIR
        rm -r $OUTPUT_DIR/data
        """
