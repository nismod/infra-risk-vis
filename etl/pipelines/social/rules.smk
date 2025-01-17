rule download_rwi_meta:
    """
    Meta Relative Wealth Index

    Description:
    https://dataforgood.facebook.com/dfg/tools/relative-wealth-index

    Catalogue page: https://data.humdata.org/dataset/relative-wealth-index


    Citation:

    Microestimates of wealth for all low- and middle-income countries Guanghua
    Chi, Han Fang, Sourav Chatterjee, Joshua E. Blumenstock Proceedings of the
    National Academy of Sciences Jan 2022, 119 (3) e2113658119; DOI:
    10.1073/pnas.2113658119

    License: CC-BY-NC 4.0
    """
    output:
        links = "raster/raw/social/meta_rwi/rwi_links.txt",
        meta = "raster/raw/social/meta_rwi/rwi.csv",
    run:
        import pandas as pd
        from hdx.data.dataset import Dataset
        from hdx.api.configuration import Configuration

        Configuration.create(hdx_site="prod", user_agent="nismod/infra-risk-vis", hdx_read_only=True)
        ds = Dataset.read_from_hdx("relative-wealth-index")
        rs = ds.get_resources()
        meta = pd.DataFrame(rs)
        meta[['id', 'name', 'description', 'created', 'download_url']].to_csv(output.meta, index=False)
        meta[['download_url']].to_csv(output.links, header=False, index=False)

rule download_rwi:
    output:
        archive = "raster/raw/social/rwi/relative-wealth-index-april-2021.zip",
        csv = "raster/raw/social/rwi/relative-wealth-index-april-2021.csv",
    shell:
        """
        output_dir=$(dirname {output.csv})

        wget -nc \
            https://data.humdata.org/dataset/76f2a2ea-ba50-40f5-b79c-db95d668b843/resource/de2f953e-940c-43bb-b1f8-4d02d28124b5/download/relative-wealth-index-april-2021.zip \
            --directory-prefix=$(dirname {output.archive})

        unzip -j -n {output.archive} -d $output_dir

        pushd $output_dir
            awk '(NR == 1) || (FNR > 1)' *_relative_wealth_index.csv > $(basename {output.csv})
        popd
        """

rule process_rwi_points:
    input:
        csv = "raster/raw/social/rwi/relative-wealth-index-april-2021.csv",
    output:
        parquet = "raster/raw/social/rwi/relative-wealth-index-april-2021.parquet",
        gpkg = "raster/raw/social/rwi/relative-wealth-index-april-2021.gpkg",
    run:
        import pandas as pd
        import geopandas as gpd
        df = pd.read_csv(input.csv)
        gdf = gpd.GeoDataFrame(
            data=df[["rwi", "error"]],
            geometry=gpd.points_from_xy(df.longitude, df.latitude),
            crs="EPSG:4326"
        )
        # reproject to 3857 (point coordinates are provided in lat/lon but the
        # underlying regular grid seems to be based on Web Mercator via Quadkeys)
        gdf.to_crs("EPSG:3857", inplace=True)
        gdf.to_parquet(output.parquet)
        gdf.to_file(output.gpkg, engine="pyogrio")

rule process_rwi_points_tiff:
    """Rasterise points to a grid data format

    -tr specifies x/y resolution, guessed from most common difference between point
    coordinates:

    def guess_resolution(coords):
        coords.sort()
        diffs = np.diff(coords)
        diffs.sort()
        # vals, counts = np.unique(diffs, return_counts=True)
        counts, vals = np.histogram(diffs)
        return vals[np.argmax(counts)]
    print("x res:", guess_resolution(gdf.geometry.x.unique())
    print("y res:", guess_resolution(gdf.geometry.y.unique())
    """
    input:
        gpkg = "raster/raw/social/rwi/relative-wealth-index-april-2021.gpkg",
    output:
        tiff = "raster/raw/social/rwi.tif",
    shell:
        """
        gdal_rasterize \
            -a rwi \
            -init -999 \
            -a_nodata -999 \
            -tr 2445.9786434 2445.96770335 \
            -ot Float64 \
            -of GTiff \
            {input.gpkg} \
            {output.tiff}
        """

rule download_worldpop_tiff:
    """
    WorldPop Age and Sex Structures in 2020, 1km Unconstrained Global Mosaic

    Link: https://hub.worldpop.org/geodata/summary?id=24798

    Estimates of total number of people per grid square broken down by sex and
    age groupings (including 0-1 and by 5-year up to 80+) in 2020

    The dataset is available to download in Geotiff format at a resolution of
    1km. The projection is Geographic Coordinate System, WGS84. The units are
    estimated number of male/female in each age group per grid square. The
    mapping approach is Pezzulo, C. et al. Sub-national mapping of population
    pyramids and dependency ratios in Africa and Asia. Sci. Data 4:170089
    doi:10.1038/sdata.2017.89 (2017)

    Filenames: Example - global_f_5_2020_1km.tif People per pixel (PPP) for
    female age group 5 to 9 years (f_05) for year 2020. For other datasets, m =
    male, 00 = age group 0 to 12months, 01 = age group 1 to 4 years, 80 = age 80
    years and over.

    License: CC-BY 4.0 International

    Data Citation:

    WorldPop (www.worldpop.org - School of Geography and Environmental Science,
    University of Southampton; Department of Geography and Geosciences,
    University of Louisville; Departement de Geographie, Universite de Namur)
    and Center for International Earth Science Information Network (CIESIN),
    Columbia University (2018). Global High Resolution Population Denominators
    Project - Funded by The Bill and Melinda Gates Foundation (OPP1134076).
    https://dx.doi.org/10.5258/SOTON/WP00654

    Method Citation:

    Pezzulo, C. et al. Sub-national mapping of population pyramids and
    dependency ratios in Africa and Asia. Sci. Data 4:170089
    doi:10.1038/sdata.2017.89 (2017)
    """
    output:
        tiff = "raster/raw/social/worldpop/global_{GENDER}_{AGE_GROUP}_{YEAR}_1km.tif"
    shell:
        """
        base_url=https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/{wildcards.YEAR}/0_Mosaicked/global_mosaic_1km
        wget -nc \
            $base_url/global_{wildcards.GENDER}_{wildcards.AGE_GROUP}_{wildcards.YEAR}_1km.tif \
            --directory-prefix=$(dirname {output.tiff})
        """

rule download_worldpop:
    input:
        tiffs = expand(
            "raster/raw/social/worldpop/global_{GENDER}_{AGE_GROUP}_{YEAR}_1km.tif",
            GENDER=['m', 'f'],
            AGE_GROUP=[0, 1] + list(range(5, 85, 5)),
            YEAR=[2020],
        ),

rule download_gridded_hdi:
    """
    Global High-Resolution Estimates of the United Nations Human Development
    Index

    Global estimates of United Nations Human Development Index (HDI) on a global
    0.1 degree grid. Developed using a generalizable machine learning
    downscaling technique based on satellite imagery that allows for training
    and prediction with observations of arbitrary shape and size. This
    downscales the national HDI, which is a multi-dimensional index used for
    measuring national development, incorporating measures of income, education
    and health.

    Citation

    Sherman, L., et al. 2023. Global High-Resolution Estimates of the United
    Nations Human Development Index Using Satellite Imagery and
    Machine-learning. Working Paper Series. 31044. National Bureau of Economic
    Research. DOI: 10.3386/w31044 Available online:
    http://www.nber.org/papers/w31044

    License: MIT

    Description URL: https://www.mosaiks.org/hdi

    Data URL: https://github.com/Global-Policy-Lab/hdi_downscaling_mosaiks
    """
    output:
        archive = "raster/raw/social/hdi/hdi_raster_predictions_V2.0.zip",
        tiff = "raster/raw/social/hdi/hdi_raster_predictions.tif"
    shell:
        """
        wget -nc \
            https://github.com/Global-Policy-Lab/hdi_downscaling_mosaiks/raw/master/data/preds/hdi_raster_predictions_V2.0.zip \
            --directory-prefix=$(dirname {output.archive})

        unzip -j -n {output.archive} $(basename {output.tiff}) -d $(dirname {output.tiff})
        """

rule clip_raster:
    """
    Clip raster extent to window defined by `raster_bounds` in config.
    """
    input:
        "raster/raw/social/{KEY}.tif"
    output:
        temp("raster/clip/social/{KEY}.tif")
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
            -t_srs EPSG:4326 \
            -te {params.bounds} \
            -of GTiff \
            {input} \
            {output}
        """
