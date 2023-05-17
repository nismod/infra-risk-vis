# multistage build of nginx image

# stage 1 of 2
FROM node:18-bullseye AS builder
WORKDIR /app

# install npm packages -- do this before copying everything else
# to enable more efficient caching by docker (we don't reinstall packages when
# we've only changed application code)
RUN apt-get update && apt-get install -y build-essential python3 python-is-python3
COPY package.json package-lock.json ./
RUN npm ci
COPY public public
COPY src src
COPY tsconfig.json tsconfig-typings.json .eslintrc.js ./
# build to static site
RUN npm run build

# stage 2 of 2
FROM nginx:alpine
WORKDIR /usr/share/nginx/html
# remove default nginx static assets
RUN rm -rf ./*
# copy in static assets from builder stage
COPY --from=builder /app/build .
# and copy in the nginx config...
# delete conf.d as it will contain default.conf which will conflict
RUN rm -rf /etc/nginx/conf.d
COPY nginx /etc/nginx
# SSL Certs mounted externally

# containers run nginx with global directives and daemon off
ENTRYPOINT ["nginx", "-g", "daemon off;"]
