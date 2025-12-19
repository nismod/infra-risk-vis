"""Core raster processing pipeline to create Cloud-Optimised GeoTIFFs

N.B. Dataset specific rules are located in pipelines/<dataset>.
"""


rule generate_terracotta_colourmap:
    input:
        csv = "pipelines/{DATASET}/colourmap.csv",
        script = "scripts/legend_to_tc_colourmap.py"
    output:
        json = "pipelines/{DATASET}/colourmap.json"
    shell:
        """
        python {input.script} {input.csv} {output.json}
        """


rule clip_raster:
    """
    Clip raster extent to window defined by `raster_bounds` in config.
    """
    input:
        "raster/raw/{DATASET}/{KEY}.tif"
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
            -t_srs EPSG:4326 \
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
        mem_mb=20000,
    priority:
        90,
    shell:
        """
        terracotta optimize-rasters \
            -o $(dirname {output}) \
            --overwrite \
            --in-memory \
            --reproject \
            --nproc 1 \
            --resampling-method nearest \
            {input}
        """
