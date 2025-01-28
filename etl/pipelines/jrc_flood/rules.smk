rule download:
    output:
        zip="raster/raw/jrc_flood/floodMapGL_rp{RP}y.zip"
    shell:
        """
        output_dir=$(dirname {output.zip})

        wget -q -nc \
            --directory-prefix=$output_dir \
            https://cidportal.jrc.ec.europa.eu/ftp/jrc-opendata/FLOODS/GlobalMaps/floodMapGL_rp{wildcards.RP}y.zip
        """

rule extract:
    input:
        tiff="raster/raw/jrc_flood/floodMapGL_rp{RP}y.zip"
    output:
        tiff="raster/raw/jrc_flood/floodMapGL_rp{RP}y.tif"
    shell:
        """
        output_dir=$(dirname {output.tiff})
        unzip $output_dir/floodMapGL_rp{wildcards.RP}y.zip floodMapGL_rp{wildcards.RP}y.tif -d $output_dir
        """

rule all:
    input:
        tiffs=expand("raster/raw/jrc_flood/floodMapGL_rp{RP}y.tif", RP=[10, 20, 50, 100, 200, 500])
