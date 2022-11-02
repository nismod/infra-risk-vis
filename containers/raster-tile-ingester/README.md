# Raster Tileserver Ingester

Utilities for managing Cloud-Optimised Tiff files in a Terracotta MySQL database.

### Environment

The following environment variables are required:

```
#### Terracotta Env
TC_DRIVER_PATH=mysql://{mysql user}:{mysql passwd}@{mysql host}
TC_DRIVER_PROVIDER=mysql
TC_PNG_COMPRESS_LEVEL=0
TC_RESAMPLING_METHOD="nearest"
TC_REPROJECTION_METHOD="nearest"

#### Gri Backend Env - for managing entries in the internal API tileserver
BACKEND_HOST=
BACKEND_PORT=
API_TOKEN=
```