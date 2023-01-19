from uuid import UUID

import pytest
from shared.infrastructure.events import InternalDomainEventDispatcher
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from substances_core.application.syncs.commands.execute import (
    ExecuteSyncCommand,
    ExecuteSyncCommandHandler,
)
from substances_core.application.syncs.commands.internal.extract import (
    ExtractSyncCommand,
    ExtractSyncCommandHandler,
)
from substances_postgresql.uow import PostgresqlUnitOfWork
from substances_rabbitmq.queue import RabbitMessageQueue, RabbitMessageQueueOptions
from substances_rest.factories import SubstanceExtractorFactory


@pytest.fixture()
def database_engine():
    engine = create_engine(
        "postgresql://localhost/substances?user=admin&password=admin"
    )
    return engine


@pytest.fixture()
def database_session():
    engine = create_engine(
        "postgresql://localhost/substances?user=admin&password=admin"
    )
    session = sessionmaker(autocommit=False, autoflush=False, bind=engine)()
    try:
        yield session
    finally:
        session.close()


@pytest.fixture()
def unit_of_work(database_session):
    return PostgresqlUnitOfWork(database_session)


@pytest.fixture()
def domain_event_dispatcher():
    return InternalDomainEventDispatcher()


@pytest.fixture()
def rabbit_message_queue():
    options = RabbitMessageQueueOptions(queue="test", host="localhost", port=5672)
    return RabbitMessageQueue(options)


@pytest.fixture()
def substance_extractor_factory():
    return SubstanceExtractorFactory()


def test_execute_command_handler(unit_of_work, rabbit_message_queue):
    # Arrange
    handler = ExecuteSyncCommandHandler(unit_of_work, rabbit_message_queue)
    command = ExecuteSyncCommand(provider="PUBCHEM")

    # Act
    result = handler.execute(command)

    # Assert
    assert result is not None


def test_checkin_command_handler(
    unit_of_work,
    rabbit_message_queue,
    substance_extractor_factory,
):
    # Arrange
    handler = ExtractSyncCommandHandler(
        unit_of_work,
        rabbit_message_queue,
        substance_extractor_factory,
    )
    command = ExtractSyncCommand(sync_id=UUID("999c1550-cb56-4a18-b09b-bba44b26c089"))

    # Act
    result = handler.execute(command)

    # Assert
    assert result is not None


def test_checkout_command_handler():
    # Arrange

    pass
