# Vector Tile Server

`tileserver-gl-light` is used to serve MBTiles containing vector data.

### Docker

We build on the
[maptiler/tileserver-gl](https://hub.docker.com/r/maptiler/tileserver-gl) image
with command customised in each docker compose file. This mounts `fonts`,
`data`, and `config.json`. The container serves vector tiles, listening on the
port and base URL specified in the `command` block.
