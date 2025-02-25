rule download_all:
    output:
        txt="raster/raw/jrc_flood/README.txt"
    shell:
        """
        output_dir=$(dirname {output.zip})

        wget --recursive --no-parent --continue --no-clobber \
            --no-host-directories \
            --cut-dirs=4 \
            --directory-prefix=$output_dir \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/CEMS-GLOFAS/flood_hazard/
        """

rule combine:
    input:
        tiffs="raster/raw/jrc_flood/floodMapGL_rp{RP}y.zip"
    output:
        tiff="raster/raw/jrc_flood/jrc_2.1.0_rp{RP}.tif"
    shell:
        """
        output_dir=$(dirname {output.tiff})
        # gdal_mosaic?
        """

rule all:
    input:
        tiffs=expand("raster/raw/jrc_flood/floodMapGL_rp{RP}y.tif", RP=[10, 20, 50, 100, 200, 500])
