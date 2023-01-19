FROM mhart/alpine-node:latest

ENV PORT 3200

WORKDIR /app
COPY . .

RUN apk update \
  && npm i -g npm@latest \
  && npm i

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
  CMD curl -f http://localhost:$PORT/version || exit 1

EXPOSE ${PORT}

CMD [ "yarn", "start" ]