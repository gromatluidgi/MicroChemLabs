from fastapi import APIRouter, Depends, status
from shared.application.commands import CommandHandler
from substances_core.application.syncs.commands.execute import (
    ExecuteSyncCommand,
    ExecuteSyncCommandResult,
)
from substances_core.application.syncs.queries.get_syncs import GetSyncsQueryResult

from .dependencies import execute_sync_command_handler

router = APIRouter()


@router.get(
    "",
    status_code=status.HTTP_200_OK,
    responses={
        status.HTTP_200_OK: {"model": GetSyncsQueryResult},
    },
)
def get_syncs() -> GetSyncsQueryResult:
    """Retrieve all substance synchronizations."""
    return GetSyncsQueryResult()


@router.post(
    "",
    responses={
        status.HTTP_200_OK: {"model": ExecuteSyncCommandResult},
    },
)
def execute_sync(
    command: ExecuteSyncCommand,
    handler: CommandHandler[ExecuteSyncCommand, ExecuteSyncCommandResult] = Depends(
        execute_sync_command_handler
    ),
) -> ExecuteSyncCommandResult:
    """Create a new substance synchronization."""
    return handler.execute(command)
