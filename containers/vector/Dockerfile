FROM maptiler/tileserver-gl:v4.11.1
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

WORKDIR /

# copy in required config
COPY ./fonts /fonts
COPY ./config.json .

# make the data folder readable
# i would have thought `USER root` is redundant here, but chmod will fail without it
USER root
RUN chmod -R u+r /data

CMD ["tileserver-gl-light", "-c", "config.json", "-p", "8080", "--verbose", "-u", "http://localhost/vector"]
