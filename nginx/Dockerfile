FROM nginx:1.19.0-alpine
ARG BIGGR_SERVER_TYPE=prod

COPY ./default-$BIGGR_SERVER_TYPE.conf /etc/nginx/conf.d/default.conf
RUN [ "$BIGGR_SERVER_TYPE" == "prod" ] && mkdir -p /etc/nginx/ssl/live/biggr.org || echo "Not in production"

# The base image symlinks access.log to /dev/stdout. That is fine for
# `docker logs`, but the log exporter needs a real file it can tail and
# re-read, so replace the symlink with one. Errors still go to stderr and
# stay visible in `docker logs`.
RUN rm -f /var/log/nginx/access.log && touch /var/log/nginx/access.log
