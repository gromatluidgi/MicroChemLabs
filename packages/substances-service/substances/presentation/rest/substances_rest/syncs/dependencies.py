from fastapi import Depends
from substances_core.application.syncs.commands.execute import ExecuteSyncCommandHandler

from ..dependencies import postgres_unit_of_work, rabbit_message_queue


# Command Handlers
def execute_sync_command_handler(
    uow=Depends(postgres_unit_of_work),
    message_queue=Depends(rabbit_message_queue),
) -> ExecuteSyncCommandHandler:
    return ExecuteSyncCommandHandler(uow, message_queue)


# Query Handlers
