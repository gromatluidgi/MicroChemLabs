FROM python:3.9.5-slim as base-image

# Setup environment variables
ENV PIP_NO_CACHE_DIR=off \
    POETRY_VERSION=1.2.2

RUN apt-get update \
    && apt-get -y install libpq-dev gcc python-rdkit librdkit1 rdkit-data libxrender1 libxtst6 libxi6

# Install poetry
RUN pip install --upgrade pip
RUN pip install "poetry==$POETRY_VERSION"
# Disabled poetry venv will install package globally
RUN poetry config virtualenvs.create false

# Transfer sources
WORKDIR /home
COPY shared-kernel ./shared-kernel
COPY science-kernel ./science-kernel

COPY substances-service/substances/core ./core
COPY substances-service/substances/database/postgresql ./database/postgresql
COPY substances-service/substances/infrastructure/rabbitmq ./infrastructure/rabbitmq
COPY substances-service/substances/infrastructure/importers ./infrastructure/importers
COPY substances-service/substances/infrastructure/etl ./infrastructure/etl
COPY substances-service/substances/presentation/rest ./presentation/rest

# Setup
WORKDIR /home/shared-kernel
RUN poetry install --all-extras --without dev

WORKDIR /home/science-kernel
RUN poetry install --all-extras --without dev

WORKDIR /home/core
RUN poetry install --without dev

WORKDIR /home/database/postgresql
RUN poetry install --all-extras --without dev

WORKDIR /home/infrastructure/etl
RUN poetry install --all-extras --without dev

WORKDIR /home/infrastructure/importers
RUN poetry install --without dev

WORKDIR /home/presentation/rest
RUN poetry install --without dev

# Clean
WORKDIR /home
RUN rm -rf /shared-kernel
RUN rm -rf /science-kernel
RUN rm -rf /core
RUN rm -rf /database/postgresql
RUN rm -rf /infrastructure/etl
RUN rm -rf /infrastructure/importers

WORKDIR /home/presentation/rest

EXPOSE 5400

CMD ["poetry", "run", "uvicorn", "substances_rest.main:api", "--host", "0.0.0.0", "--port", "5400"]