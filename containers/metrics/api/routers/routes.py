# Retrieval of all country metrics data
API_ROUTE_BASE = "/metrics"

# Retrieval of GDL data and associated geojson
GDL_BASE_ROUTE = API_ROUTE_BASE + "/gdl"

# GDL geojson
GDL_BOUNDARY_BASE_ROUTE = GDL_BASE_ROUTE + "/geojson"
GDL_SUBNATIONAL_FROM_ISO_ROUTE = GDL_BOUNDARY_BASE_ROUTE + "/subnational/iso/{iso_code}"
GDL_NATIONAL_FROM_ISO_ROUTE = GDL_BOUNDARY_BASE_ROUTE + "/national/iso/{iso_code}"

# GDL countries and regions metadata
GDL_META_BASE_ROUTE = GDL_BASE_ROUTE + "/meta"
GDL_COUNTRY_META_ROUTE = GDL_META_BASE_ROUTE + "/countries"
GDL_REGION_META_ROUTE = GDL_META_BASE_ROUTE + "/regions"

# GDL annual metrics
GDL_DATA_BASE_ROUTE = GDL_BASE_ROUTE + "/data"
GDL_DATA_ISO_ROUTE = GDL_DATA_BASE_ROUTE + "/{dataset_key}" + "/{iso_code}"
GDL_DATA_EXTENT_ROUTE = GDL_DATA_BASE_ROUTE + "/{dataset_key}" + "/extent"
