"""DEM processing

Overrides generic raster processing rules to handle large DEM TIFFs.
"""

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
