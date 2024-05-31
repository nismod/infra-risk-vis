"""
These raster processing rules are reused across many datasets. They form the
core of the raster processing pipeline.

N.B. Dataset specific rules are located in pipelines/<dataset>.
"""

from pipelines.helpers import gdalwarp_bounds


rule set_zero_to_no_data:
    input:
        "raster/raw/{DATASET}/{KEY}.tif"
    output:
        temp("raster/no_data/{DATASET}/{KEY}.tif")
    resources:
        disk_mb=3000,
        mem_mb=10000,
    priority:
        70,
    shell:
        """
        NODATA=$(gdalinfo "{input}" -json | jq .bands[0].noDataValue)

        # handle case of NODATA == nan - the JSON output of gdalinfo will change
        # nan to "NaN" so we need to reverse that for gdal_calc.py
        if [ "$NODATA" == '"NaN"' ]
        then
          NODATA=nan
        fi

        if [ "$NODATA" == 'null' ]
        then
          NODATA=nan
        fi

        # replace zeros with NoData value
        gdal_calc.py \
          --quiet \
          -A "{input}" \
          --co="COMPRESS=LZW" \
          --outfile="{output}" \
          --overwrite \
          --calc="numpy.where(A<=0,$NODATA,A)" \
          --NoDataValue=$NODATA \
          --hideNoData
        """


rule clip_raster:
    """
    Clip raster extent to window defined by `raster_bounds` in config.
    """
    input:
        "raster/no_data/{DATASET}/{KEY}.tif"
    output:
        temp("raster/clip/{DATASET}/{KEY}.tif")
    params:
        bounds = config["raster_bounds"]
    resources:
        disk_mb=3000,
        mem_mb=10000,
    priority:
        80,
    shell:
        """
        gdalwarp \
            -co "COMPRESS=LZW" \
            -te {params.bounds} \
            -of GTiff \
            {input} \
            {output}
        """


rule cloud_optimise_raster:
    """
    Use terracotta to cloud optimise rasters.
    """
    input:
        "raster/clip/{DATASET}/{KEY}.tif"
    output:
        "raster/cog/{DATASET}/{KEY}.tif",
    resources:
        disk_mb=100,
        mem_mb=2000,
    priority:
        90,
    shell:
        """
        terracotta optimize-rasters \
            -o $(dirname {output}) \
            --overwrite \
            --reproject \
            --nproc -1 \
            --resampling-method nearest \
            {input}
        """
