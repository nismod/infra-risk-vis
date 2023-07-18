rule download_motorised:
    """
    Download motorised travel time to nearest healthcare facility.
    """
    output:
        archive = temp("raster/raw/traveltime_to_healthcare/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.zip"),
        raster = "raster/raw/traveltime_to_healthcare/202001_Global_Motorized_Travel_Time_to_Healthcare_2019.tif"
    shell:
        """
        curl \
            'https://data.malariaatlas.org/geoserver/Accessibility/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Accessibility:202001_Global_Motorized_Travel_Time_to_Healthcare' \
            -H 'Accept-Encoding: gzip, deflate, br' \
            -H 'Connection: keep-alive' \
            --output {output.archive}

        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """


rule download_walking:
    """
    Download walking travel time to nearest healthcare facility.
    """
    output:
        archive = temp("raster/raw/traveltime_to_healthcare/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.zip"),
        raster = "raster/raw/traveltime_to_healthcare/202001_Global_Walking_Only_Travel_Time_To_Healthcare_2019.tif"
    shell:
        """
        curl \
            'https://data.malariaatlas.org/geoserver/Accessibility/ows?service=CSW&version=2.0.1&request=DirectDownload&ResourceId=Accessibility:202001_Global_Walking_Only_Travel_Time_To_Healthcare' \
            -H 'Accept-Encoding: gzip, deflate, br' \
            -H 'Connection: keep-alive' \
            --output {output.archive}

        unzip {output.archive} $(basename {output.raster}) -d $(dirname {output.raster})
        """