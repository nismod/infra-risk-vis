rule download_gdl:
    output:
        xlsx="vector/raw/gdl/GDL Codes V6.5.xlsx",
        zip="vector/raw/gdl/GDL Shapefiles V6.5.zip",
        csv="vector/raw/gdl/Subnational HDI Data v8.3.csv",
        meta_csv="vector/raw/gdl/SHDI-SGDI-8.3-Vardescription.csv",

        # from previous release:
        json_large="vector/raw/gdl/gdl_6.4_large_vis_0.1.json",
        json_small="vector/raw/gdl/gdl_v6.4_national_small.json",
        csv_old="vector/raw/gdl/shdi_sgdi_total_8.0.csv",
    shell:
        """
        pushd $(dirname '{output.zip}')
            zenodo_get -w links_17467221.txt --record=17467221
            wget -nc -i links_17467221.txt
            md5sum -c md5sums.txt

            # previous release:
            zenodo_get -w links_14868935.txt --record=14868935
            wget -nc -i links_14868935.txt
            md5sum -c md5sums.txt
        popd
        """
