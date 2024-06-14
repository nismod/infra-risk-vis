rule download:
    """
    Download files from remote location.
    """
    output:
        zip="raster/raw/storm/STORM_FIXED_RETURN_PERIODS.zip",
    shell:
        """
        OUTPUT_DIR=$(dirname {output.zip})

        pushd $OUTPUT_DIR
            zenodo_get -w links.txt --record=10931452
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """

rule unpack:
    input:
        zip=rules.download.output.zip,
    output:
        tiff="raster/raw/storm/{KEY}.tif",
    shell:
        """
        unzip -n {input.zip} $(basename {output.tiff}) -d $(dirname {output.tiff})
        """
