from fastapi import Depends
from substances_core.application.substances.commands.import_string import (
    ImportSubstanceCommandHandler,
)
from substances_core.application.syncs.commands.execute import ExecuteSyncCommandHandler
from substances_importers.factory import DefaultSubstanceImporterFactory

from ..dependencies import postgres_unit_of_work, rabbit_message_queue


def substance_importer_factory():
    return DefaultSubstanceImporterFactory()


# Command Handlers
def execute_sync_command_handler(
    uow=Depends(postgres_unit_of_work),
    message_queue=Depends(rabbit_message_queue),
) -> ExecuteSyncCommandHandler:
    return ExecuteSyncCommandHandler(uow, message_queue)


def import_substance_from_string_command_handler(
    importer=Depends(substance_importer_factory), uow=Depends(postgres_unit_of_work)
) -> ImportSubstanceCommandHandler:
    return ImportSubstanceCommandHandler(uow, importer)


# Query Handlers
