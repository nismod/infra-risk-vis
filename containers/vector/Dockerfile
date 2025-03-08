FROM maptiler/tileserver-gl:v4.12.0
LABEL maintainer="frederick.thomas@ouce.ox.ac.uk"

# used in run.sh - override through docker compose if needed
# this URL will be contained in the responses from the tileserver
# it should be a URL accessible from the client
ENV TILESERVER_PUBLIC_URL="http://localhost/vector"

WORKDIR /

# copy in required config
COPY ./fonts /fonts
COPY ./config.json .

# make the data folder readable
# i would have thought `USER root` is redundant here, but chmod will fail without it
USER root
RUN chmod -R u+r /data

CMD ["./run.sh"]
