import logging
import threading

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from .api import api_router
from .config import settings
from .graphql.router import graphql_router
from .grpc.server import start_grpc_server
from .logging import configure_logging
from .queues.consumer import RabbitConsumer

log = logging.getLogger(__name__)

configure_logging()

api = FastAPI(
    title="Subtances Microservice API",
    description="Subtances Microservice project",
    version="1.0.0",
)


origins = ["http://localhost", "http://localhost:8080"]
api.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

api.include_router(api_router, prefix="/api/v1")
api.include_router(graphql_router)  # TODO: remove and use standalone server


@api.on_event("startup")
def startup_rabbit_consumers():
    for _ in range(settings.RABBITMQ_CONCURENT_CONSUMER):
        RabbitConsumer.run()


# TODO: remove and use standalone server
@api.on_event("startup")
def startup_grpc():
    threading.Thread(target=start_grpc_server, args=[]).start()
