# Traefik Reverse Proxy

[Traefik](https://doc.traefik.io/traefik/user-guides/docker-compose/basic-example/)
is used in both dev and prod setups as a reverse proxy and TLS certificate
managment helper.

## Development Setup

Localhost certificates and TLS configuration are included here.

## Production Setup

LetsEncrypt is configured as production certificate manager, with certificates
to be mounted in the production volume.
