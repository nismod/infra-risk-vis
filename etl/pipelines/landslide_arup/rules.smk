
rule download_metadata:
    output:
        json="raster/raw/landslide_arup/metadata.json",
    shell:
        """
        wget --output-document={output.json} \
            https://datacatalogapi.worldbank.org/ddhxext/DatasetDownload?dataset_unique_id=0037584&version_id=
        """

rule download_landslides:
    input:
        json="raster/raw/landslide_arup/metadata.json",
    output:
        "raster/raw/landslide_arup/global-landslide-hazard-map-report.pdf",
        "raster/raw/landslide_arup/ls_eq_tiled.tif",
        "raster/raw/landslide_arup/LS_RF_Mean_1980-2018_COG.tif",
        "raster/raw/landslide_arup/LS_RF_Median_1980-2018_COG.tif",
        "raster/raw/landslide_arup/LS_TH_COG.tif",
    run:
        import json
        import os
        import zipfile
        from pathlib import Path
        import requests

        with open(input.json, 'r') as fh:
            meta = json.load(fh)

        for file_meta in meta["resources"]:
            fname = file_meta["distribution"]["file_name"]
            url = file_meta["distribution"]["url"]
            out_dir = "raster/raw/landslide"
            out_file = os.path.join(out_dir, fname)

            if Path(out_file).exists():
                print("Skipped downloading", fname)
            else:
                print("Downloading", url)
                r = requests.get(url)
                with open(out_file, 'wb') as fd:
                    for chunk in r.iter_content(chunk_size=1024):
                        fd.write(chunk)

            if ".zip" in fname:
                print("Extracting zip", out_file)
                with zipfile.ZipFile(out_file, 'r') as zh:
                    zh.extractall(out_dir)


rule ingest_categorical_raster:
    """
    Custom ingestion rule for susceptibility, LS_TH_COG.tif as raster is
    categorical and needs a mapping from raster integer value to colour and
    classification.
    """
    input:
        raster = "raster/cog/landslide_arup/LS_TH_COG.tif",
        legend = "pipelines/landslide_arup/colourmap.csv",
        metadata = "pipelines/landslide_arup/metadata.json",
    params:
        key_values = ["susceptibility"]
    output:
        flag = "raster/ingest/landslide_arup.categorical.flag"
    script:
        "../../scripts/ingest_categorical.py"
