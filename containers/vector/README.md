# Vector Tile Server

`tileserver-gl-light` is used to serve MBTiles containing vector-data.

### Docker

`Dockerfile-dev` is used for a locally hosted setup (see `docker-compose-dev.yaml`)

`Dockerfile-prod` is used for a deployed setup

The only difference between these Dockerfiles is the `-u` (`--public-url`) flag, which ensures the tiles and configuration json is served with the correct path (i.e. behind the reverse-proxy)