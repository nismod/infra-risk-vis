# multistage build of nginx image

# stage 1 of 2

# name the node stage "builder"
FROM node:16 AS builder
# set working directory
WORKDIR /app
# copy all files from current directory to working dir in image
COPY . .
# install node modules and build assets
RUN npm install
RUN npm run build

# stage 2 of 2

FROM nginx:alpine
# Set working directory to nginx asset directory
WORKDIR /usr/share/nginx/html
# remove default nginx static assets
RUN rm -rf ./*
# copy static assets from builder stage
COPY --from=builder /app/build .
COPY jamaica.infrastructureresilience.org /etc/nginx/sites-available/

# Containers run nginx with global directives and daemon off
ENTRYPOINT ["nginx", "-g", "daemon off;"]
