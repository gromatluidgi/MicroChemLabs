"""Lambda Dependencies"""
import logging
import os

from substances_core.application.syncs.commands.internal.transform import (
    TransformSyncCommandHandler,
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

from .factories import SubstanceTransformerFactory

LOGGER = logging.getLogger(__name__)


def substance_transformer_factory():
    return SubstanceTransformerFactory()


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


# Command Handlers
def transform_sync_command_handler(
    uow=dynamo_unit_of_work(), transformer_factory=substance_transformer_factory()
) -> TransformSyncCommandHandler:
    return TransformSyncCommandHandler(uow, transformer_factory)
