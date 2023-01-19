"""Lambda Dependencies"""
import logging
import os

from substances_core.application.syncs.commands.internal.load import (
    LoadSyncCommandHandler,
)
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_dynamodb.uow import DynamoUnitOfWork

from .factories import SubstanceLoaderFactory

LOGGER = logging.getLogger(__name__)


def substance_loader_factory():
    return SubstanceLoaderFactory()


# Transaction Management
def dynamo_unit_of_work():
    table_names = {
        SyncReadRepository: os.environ["SYNCS_TABLE"],
        SyncWriteRepository: os.environ["SYNCS_TABLE"],
        SubstanceReadRepository: os.environ["SUBSTANCES_TABLE"],
        SubstanceWriteRepository: os.environ["SUBSTANCES_TABLE"],
    }
    uow = DynamoUnitOfWork(table_names)
    return uow


# Command Handler
def load_sync_command_handler(
    uow=dynamo_unit_of_work(), factory=substance_loader_factory()
) -> LoadSyncCommandHandler:
    return LoadSyncCommandHandler(uow, factory)
