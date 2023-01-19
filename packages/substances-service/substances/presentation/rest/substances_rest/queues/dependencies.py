"""Common dependencies"""
import logging
from functools import lru_cache

from shared.infrastructure.commands import InternalCommandDispatcher
from shared.infrastructure.events import InternalDomainEventDispatcher
from shared.infrastructure.queue import MessageQueue
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, scoped_session, sessionmaker
from substances_core.application.syncs.commands.internal.extract import (
    ExtractSyncCommand,
    ExtractSyncCommandHandler,
)
from substances_core.application.syncs.commands.internal.load import (
    LoadSyncCommand,
    LoadSyncCommandHandler,
)
from substances_core.application.syncs.commands.internal.transform import (
    TransformSyncCommand,
    TransformSyncCommandHandler,
)
from substances_core.application.syncs.events.checkout import SyncCheckoutEventHandler
from substances_core.domain.syncs.events import SyncCheckedOutEvent
from substances_postgresql.uow import PostgresqlUnitOfWork
from substances_rabbitmq.queue import RabbitMessageQueue, RabbitMessageQueueOptions

from ..config import settings
from ..factories import (
    SubstanceExtractorFactory,
    SubstanceLoaderFactory,
    SubstanceTransformerFactory,
)

LOGGER = logging.getLogger(__name__)


@lru_cache()
def get_database_engine():
    LOGGER.info("Create Queue model database engine.")
    return create_engine(str(settings.DB_URI))


def postgres_context(engine=get_database_engine()):
    session: Session = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=engine)
    )
    try:
        yield session
    finally:
        session.close()


# ETL (Extract - Load - Transform) Processors


def substance_extractor_factory():
    return SubstanceExtractorFactory()


def substance_transformer_factory():
    return SubstanceTransformerFactory()


def substance_loader_factory():
    return SubstanceLoaderFactory()


# Messaging
def domain_event_dispatcher(
    session_factory=postgres_context,
) -> InternalDomainEventDispatcher:
    dispatcher = InternalDomainEventDispatcher()
    dispatcher.register_handler(
        SyncCheckedOutEvent,
        SyncCheckoutEventHandler(
            PostgresqlUnitOfWork(session_factory, dispatcher),
            substance_loader_factory(),
        ),
    )
    return dispatcher


# Transaction Management
def postgres_unit_of_work(
    session=postgres_context(), dispatcher=domain_event_dispatcher()
):
    return PostgresqlUnitOfWork(next(session), dispatcher)


def rabbit_message_queue() -> MessageQueue:
    message_queue = RabbitMessageQueue(
        RabbitMessageQueueOptions(
            host=settings.RABBITMQ_HOST, port=settings.RABBITMQ_PORT, queue="substances"
        )
    )
    return message_queue


def internal_command_dispatcher(
    uow=postgres_unit_of_work(),
) -> InternalCommandDispatcher:
    dispatcher = InternalCommandDispatcher()
    dispatcher.register_handler(
        ExtractSyncCommand,
        ExtractSyncCommandHandler(
            uow, substance_extractor_factory(), rabbit_message_queue()
        ),
    )
    dispatcher.register_handler(
        TransformSyncCommand,
        TransformSyncCommandHandler(
            uow, substance_transformer_factory(), rabbit_message_queue()
        ),
    )
    dispatcher.register_handler(
        LoadSyncCommand,
        LoadSyncCommandHandler(uow, substance_loader_factory()),
    )
    return dispatcher
