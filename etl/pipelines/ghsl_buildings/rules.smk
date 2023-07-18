from pipelines.helpers import gdalwarp_bounds

configfile: "../../config.yml"


rule download_all_building_types:
    """
    Download JRC built up area layer.
    """
    output:
        archive = temp("raster/raw/ghsl_buildings/all_types_E{EPOCH}_R{RELEASE}_archive.zip"),
        raster = "raster/raw/ghsl_buildings/GHS_BUILT_S_E{EPOCH}_GLOBE_R{RELEASE}_54009_1000_V1_0.tif"
    shell:
        """
        wget \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R{wildcards.RELEASE}/GHS_BUILT_S_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000/V1-0/GHS_BUILT_S_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000_V1_0.zip \
            --output-document={output.archive}
        
        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule download_non_residential_type:
    """
    Download JRC non-residential built up area layer.
    """
    output:
        archive = temp("raster/raw/ghsl_buildings/non_residential_E{EPOCH}_R{RELEASE}_archive.zip"),
        raster = "raster/raw/ghsl_buildings/GHS_BUILT_S_NRES_E{EPOCH}_GLOBE_R{RELEASE}_54009_1000_V1_0.tif"
    shell:
        """
        wget \
            https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_BUILT_S_GLOBE_R{wildcards.RELEASE}/GHS_BUILT_S_NRES_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000/V1-0/GHS_BUILT_S_NRES_E{wildcards.EPOCH}_GLOBE_R{wildcards.RELEASE}_54009_1000_V1_0.zip \
            --output-document={output.archive}
        
        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule clip_and_reproject_raster:
    """
    Reproject from Mollweide to WGS84 and clip to bounds.
    """
    input:
        raster = "raster/no_data/ghsl_buildings/{KEY}.tif",
    params:
        bounds = gdalwarp_bounds(config["raster_bounds"])
    output:
        "raster/clip/ghsl_buildings/{KEY}.tif"
    resources:
        mem_mb=10000
    priority:
        80,
    shell:
        """
        gdalwarp \
            -t_srs EPSG:4326 \
            -of GTiff \
            -co COMPRESS=LZW \
            -te {params.bounds} \
            {input.raster} \
            {output}
        """