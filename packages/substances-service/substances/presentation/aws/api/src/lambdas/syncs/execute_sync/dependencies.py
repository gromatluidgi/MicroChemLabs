"""Lambda Dependencies"""
import logging
import os

from shared.infrastructure.queue import MessageQueue
from substances_aws_queue.queue import SQSMessageQueue
from substances_core.application.syncs.commands.execute import ExecuteSyncCommandHandler
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_dynamodb.uow import DynamoUnitOfWork

LOGGER = logging.getLogger(__name__)


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
def execute_sync_command_handler(
    uow: DynamoUnitOfWork = dynamo_unit_of_work(),
) -> ExecuteSyncCommandHandler:
    return ExecuteSyncCommandHandler(uow)
