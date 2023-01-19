from shared.domain.specification import Specification
from substances_core.domain.susbtances.repository import SubstanceReadRepository
from substances_core.domain.susbtances.substance import Substance


class UniqueSubstanceSpecification(Specification[Substance]):
    def __init__(self, repository: SubstanceReadRepository) -> None:
        self._repository = repository

    def is_satisfied_by(self, input: Substance) -> bool:
        substance = self._repository.find_by_inchi_key(input.inchi_key)
        return substance is None
