FROM python:3.9.5-slim as base-image

ENV PORT=50051 \
    POETRY_VERSION=1.2.2

RUN apt-get update \
    && apt-get -y install libpq-dev gcc python3-pip

# Install poetry
RUN pip install --upgrade pip
RUN pip install "poetry==$POETRY_VERSION"
# Disabled poetry venv will install package globally
RUN poetry config virtualenvs.create false

WORKDIR /app
COPY ./substances-service/substances/presentation/grpc .
RUN poetry install


EXPOSE ${PORT}

CMD ["python", "server/main.py"]