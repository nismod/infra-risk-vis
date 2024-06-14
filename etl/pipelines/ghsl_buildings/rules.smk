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
