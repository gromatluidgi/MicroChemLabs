import logging
from typing import Any, List

from sciences.chemistry.molecules import MoleculeUtils
from substances_core.domain.susbtances.objects import CanonicalSmiles, Formula, InchiKey
from substances_core.domain.susbtances.services.importer import (
    ImporterDataType,
    SubstanceImporter,
)
from substances_core.domain.susbtances.substance import Substance
from substances_importers.resolvers.cactus_resolver import CactusResolver

LOGGER = logging.getLogger(__name__)


class SmilesSubtanceImporter(SubstanceImporter):
    def __init__(
        self,
        source_type: ImporterDataType,
    ) -> None:
        self._source_type = source_type

    def load(self, data: Any) -> List[Substance]:
        if self._source_type == ImporterDataType.FILE:
            pass
        elif self._source_type == ImporterDataType.STRING:
            # TODO: guard
            return self._parse_smiles_from_string(data)
        elif self._source_type == ImporterDataType.BINARY:
            pass

        raise RuntimeError

    def _parse_smiles_from_string(self, chain: str) -> List[Substance]:
        substances = []
        lines = chain.split("\n")
        for line in lines:
            try:
                canonical_smiles = CanonicalSmiles(
                    MoleculeUtils.get_canonical_smiles(line)
                )
                molecule = MoleculeUtils.get_mol_from_smiles(str(canonical_smiles))
                inchi_key = InchiKey(
                    MoleculeUtils.get_inchi_key_from_smiles(str(canonical_smiles))
                )
                iupac_name = CactusResolver.resolve_iupac_name(str(canonical_smiles))
                formula = Formula(MoleculeUtils.get_molecular_formula(molecule))
                mol_weight = MoleculeUtils.get_mol_weight(molecule)

                substance = Substance(
                    inchi_key=inchi_key,
                    formula=formula,
                    mol_weight=mol_weight,
                    canonical_smiles=canonical_smiles,
                    iupac_name=iupac_name,
                )

                substances.append(substance)
            except Exception as err:
                LOGGER.error(err)

        return substances
