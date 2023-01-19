"""Lambda Dependencies"""
import logging
import os

from shared.infrastructure.queue import MessageQueue
from substances_aws_queue.queue import SQSMessageQueue
from substances_core.application.syncs.commands.internal.extract import (
    ExtractSyncCommandHandler,
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

from .factories import SubstanceExtractorFactory  # type: ignore

LOGGER = logging.getLogger(__name__)


def substance_extractor_factory() -> SubstanceExtractorFactory:
    return SubstanceExtractorFactory()


# Transaction Management
def dynamo_unit_of_work() -> DynamoUnitOfWork:
    table_names = {
        SyncReadRepository: os.environ["SYNCS_TABLE"],
        SyncWriteRepository: os.environ["SYNCS_TABLE"],
        SubstanceReadRepository: os.environ["SUBSTANCES_TABLE"],
        SubstanceWriteRepository: os.environ["SUBSTANCES_TABLE"],
    }
    uow = DynamoUnitOfWork(table_names)
    return uow


def sqs_message_queue() -> MessageQueue:
    message_queue = SQSMessageQueue("QUEUE_NAME")
    return message_queue


# Command Handlers
def extract_sync_command_handler(
    uow: DynamoUnitOfWork = dynamo_unit_of_work(),
    extract_factory: SubstanceExtractorFactory = substance_extractor_factory(),
) -> ExtractSyncCommandHandler:
    return ExtractSyncCommandHandler(uow, extract_factory)
