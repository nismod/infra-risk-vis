rule download_1000m_organic_carbon_stocks:
    """
    Download organic carbon stock raster at 1000m resolution and reproject to WGS84.
    """
    output:
        homolosine = temp("raster/raw/exposure_nature/ocs_0-30cm_mean_1000_homolosine.tif"),
        WGS84 = "raster/raw/exposure_nature/ocs_0-30cm_mean_1000.tif"
    shell:
        """
        wget \
        https://files.isric.org/soilgrids/latest/data_aggregated/1000m/ocs/ocs_0-30cm_mean_1000.tif \
        --tries 3 \
        --output-document={output.homolosine}

        gdalwarp \
            -of Gtiff \
            -co COMPRESS=LZW \
            -t_srs EPSG:4326 \
            {output.homolosine} \
            {output.WGS84}
        """


rule download_3arcsec_biodiversity_intactness:
    """
    Download biodiversity intactness data. Unzip and translate from ASCII to binary TIFF.

    N.B. The download link was generated by following "Download" and entering an email on the following page:
    https://data.nhm.ac.uk/dataset/global-map-of-the-biodiversity-intactness-index-from-newbold-et-al-2016-science/resource/8531b4dc-bd44-4586-8216-47b3b8d60e85

    The download link may expire in time.
    """
    output:
        archive = temp("raster/raw/exposure_nature/lbii.zip"),
        ascii_text = temp("raster/raw/exposure_nature/lbii.asc"),
        tiff = "raster/raw/exposure_nature/lbii.tif",
    shell:
        """
        wget \
            https://data.nhm.ac.uk/dataset/17179b71-c5e1-435c-b3db-ebc7a65d980a/resource/8531b4dc-bd44-4586-8216-47b3b8d60e85/download/lbii.zip \
            --tries 3 \
            --output-document={output.archive}

        unzip {output.archive} -d $(dirname {output.ascii_text})

        gdal_translate \
            {output.ascii_text} \
            {output.tiff} \
            -co "COMPRESS=LZW" \
            -a_srs "EPSG:4326"
        """


rule download_300m_forest_integrity_index:
    """
    Download 300m resolution forest integrity index data from Google Drive and
    rescale from [0, 10000] integer -> [0, 10] float.

    Google Drive folder is located here:
    https://drive.google.com/drive/folders/180DXlbF4dwCYhBW025YkZbhNHnrW5NlW
    """
    output:
        raw_integer = temp("raster/raw/exposure_nature/flii_integer.tif"),
        rescaled = "raster/raw/exposure_nature/flii_earth.tif"
    shell:
        """
        gdown --output {output.raw_integer} 1Bd3LxqPTSMuFRb-24Z7UMWCiPnhlTvg_ 

        gdal_translate \
            -scale 0 10000 0 10 \
            -ot Float32 \
            {output.raw_integer} \
            {output.rescaled}
        """