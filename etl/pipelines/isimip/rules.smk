rule download:
    """
    Fetch our preprocessed data from zenodo.
    """
    output:
        zip="raster/raw/isimip/lange2020_expected_occurrence.zip",
    shell:
        """
        pushd $(dirname {output.zip})
            zenodo_get -w links.txt --record=8147088
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """

rule unpack:
    input:
        zip=rules.download.output.zip,
    output:
        tiff="raster/raw/isimip/{KEY}.tif",
    shell:
        """
        unzip -n {input.zip} $(basename {output.tiff}) -d $(dirname {output.tiff})
        """
