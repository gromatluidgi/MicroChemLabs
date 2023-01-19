from typing import List

from shared.application.commands import CommandBase, CommandHandler, CommandResultBase
from shared.domain.uow import UnitOfWork
from substances_core.domain.susbtances.repository import SubstanceWriteRepository
from substances_core.domain.susbtances.services.importer import (
    ImporterDataType,
    SubstanceImporterFactory,
)
from substances_core.domain.susbtances.substance import Substance


class ImportSubstanceCommand(CommandBase):
    data: str
    format: str
    data_type: ImporterDataType = ImporterDataType.STRING


class ImportSubstanceResult(CommandResultBase):
    substances: List[Substance]


class ImportSubstanceCommandHandler(
    CommandHandler[ImportSubstanceCommand, ImportSubstanceResult]
):
    def __init__(
        self,
        unit_of_work: UnitOfWork,
        substance_importer_factory: SubstanceImporterFactory,
    ) -> None:
        self._unit_of_work = unit_of_work
        self._substance_importer_factory = substance_importer_factory

    def execute(self, command: ImportSubstanceCommand) -> ImportSubstanceResult:

        importer = self._substance_importer_factory.importer_by_format(
            ImporterDataType.STRING, command.format
        )

        substances = []

        with self._unit_of_work as uow:
            substance_write_repo: SubstanceWriteRepository = uow.get_repository(
                SubstanceWriteRepository
            )  # type: ignore

            substances = importer.load(command.data)
            substance_write_repo.batch_insert_or_update(substances)
            uow.commit()

        return ImportSubstanceResult(substances=substances)
