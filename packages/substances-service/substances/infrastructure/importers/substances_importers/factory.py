from substances_core.domain.susbtances.services.importer import (
    ImporterDataType,
    SubstanceImporter,
    SubstanceImporterFactory,
)
from substances_importers.smiles_importer import SmilesSubtanceImporter


class DefaultSubstanceImporterFactory(SubstanceImporterFactory):
    def importer_by_format(
        self, source_type: ImporterDataType, format: str
    ) -> SubstanceImporter:
        if format == "SMILES":
            return SmilesSubtanceImporter(source_type)
