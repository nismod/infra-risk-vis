"""DEM processing

Overrides generic raster processing rules to handle large DEM TIFFs.
"""

from pipelines.helpers import gdalwarp_bounds

rule download:
    output:
        "raster/raw/dem/dtm_elevation_merit.dem_m_250m_s0..0cm_2017_v1.0.tif",
        "raster/raw/dem/dtm_slope_merit.dem_m_250m_s0..0cm_2017_v1.0.tif"
    shell:
        """
        cd raster/raw/dem/

        wget -nc https://zenodo.org/record/1447210/files/dtm_elevation_merit.dem_m_250m_s0..0cm_2017_v1.0.tif
        wget -nc https://zenodo.org/record/1447210/files/dtm_slope_merit.dem_m_250m_s0..0cm_2017_v1.0.tif

        cat << EOF > md5sums.txt
e95337d45c8039141ad55ee35deaeb86  dtm_elevation_merit.dem_m_250m_s0..0cm_2017_v1.0.tif
c4fb30939620918e89de47437b673435  dtm_slope_merit.dem_m_250m_s0..0cm_2017_v1.0.tif
EOF
        md5sum -c md5sums.txt
        """


rule set_zero_to_no_data:
    input:
        "raster/raw/dem/{KEY}.tif"
    output:
        temp("raster/no_data/dem/{KEY}.tif")
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
        "raster/no_data/dem/{KEY}.tif"
    output:
        temp("raster/clip/dem/{KEY}.tif")
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
