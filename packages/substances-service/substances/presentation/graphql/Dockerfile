FROM node:16

ENV PORT 4000

WORKDIR /app
COPY ./substances-service/substances/presentation/graphql .

RUN yarn install

HEALTHCHECK --interval=30s --timeout=30s --start-period=5s --retries=3 \
  CMD curl -f http://localhost:$PORT/version || exit 1

EXPOSE ${PORT}

CMD [ "yarn", "start" ]