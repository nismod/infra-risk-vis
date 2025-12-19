"""
GEM Global Seismic Hazard map
"""

rule download:
    output:
        zip="raster/raw/gem_earthquake/GEM-GSHM_PGA-475y-rock_v2023.zip"
    shell:
        """
        pushd $(dirname {output.zip})
            zenodo_get -w links.txt --record=8409647
            wget -nc -i links.txt
            md5sum -c md5sums.txt
        popd
        """

rule unpack:
    input:
        zip=rules.download.output.zip,
    output:
        tiff="raster/raw/gem_earthquake/{KEY}.tif",
    shell:
        """
        unzip -j -n {input.zip} $(basename {output.tiff}) -d $(dirname {output.tiff})
        """
