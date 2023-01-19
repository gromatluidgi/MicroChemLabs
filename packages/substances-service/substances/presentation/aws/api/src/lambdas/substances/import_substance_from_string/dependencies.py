"""Lambda Dependencies"""
import logging

from substances_core.application.substances.commands.import_string import (
    ImportSubstanceCommandHandler,
)
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_dynamodb.uow import DynamoUnitOfWork
from substances_importers.factory import DefaultSubstanceImporterFactory

LOGGER = logging.getLogger(__name__)


def substance_importer_factory():
    return DefaultSubstanceImporterFactory()


# Transaction Management
def dynamo_unit_of_work():
    table_names = {
        SubstanceReadRepository: "substances-service-dev-substances",
        SubstanceWriteRepository: "substances-service-dev-substances",
    }
    uow = DynamoUnitOfWork(table_names)
    return uow


# Command Handlers
def import_substance_from_string_command_handler(
    uow=dynamo_unit_of_work(),
    importer=substance_importer_factory(),
) -> ImportSubstanceCommandHandler:
    return ImportSubstanceCommandHandler(uow, importer)
