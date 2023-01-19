"""Common dependencies"""
# Data Access
import logging
from functools import lru_cache

from fastapi import Depends
from shared.infrastructure.events import InternalDomainEventDispatcher
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, scoped_session, sessionmaker
from substances_postgresql.uow import PostgresqlUnitOfWork
from substances_rabbitmq.queue import RabbitMessageQueue, RabbitMessageQueueOptions

from .config import settings

LOGGER = logging.getLogger(__name__)


@lru_cache()
def get_database_engine():
    LOGGER.info("Create database engine.")
    return create_engine(str(settings.DB_URI))


def postgres_context(engine=Depends(get_database_engine)):
    session: Session = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=engine)
    )
    try:
        yield session
    finally:
        session.close()


# Messaging
def rabbit_message_queue():
    options = RabbitMessageQueueOptions(
        host=settings.RABBITMQ_HOST,
        port=settings.RABBITMQ_PORT,
        queue=settings.RABBITMQ_NAME,
    )
    return RabbitMessageQueue(options)


def domain_event_dispatcher() -> InternalDomainEventDispatcher:
    dispatcher = InternalDomainEventDispatcher()
    return dispatcher


# Transaction Management
def postgres_unit_of_work(
    session_factory=Depends(postgres_context),
    dispatcher=Depends(domain_event_dispatcher),
):
    return PostgresqlUnitOfWork(session_factory, dispatcher)
