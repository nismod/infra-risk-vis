"""
Overrides of raster processing rules to handle large DEM TIFFs.

The raw files cannot be programmatically fetched, they must be
placed in raster/raw/dem_elevation/ manually prior to running snakemake.
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
          --co="BIGTIFF=YES" \
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
        bounds = gdalwarp_bounds(config["raster_bounds"])
    resources:
        disk_mb=3000,
        mem_mb=10000,
    priority:
        80,
    shell:
        """
        gdalwarp \
            -co "COMPRESS=LZW" \
            -co "BIGTIFF=YES" \
            -te {params.bounds} \
            -of GTiff \
            {input} \
            {output}
        """
