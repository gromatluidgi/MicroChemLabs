from shared.application.queries import QueryBase, QueryHandler, QueryResultBase
from shared.domain.uow import UnitOfWork
from substances_core.domain.susbtances.objects import InchiKey
from substances_core.domain.susbtances.repository import SubstanceReadRepository
from substances_core.domain.susbtances.substance import Substance


class GetSubstanceByIdentifierQuery(QueryBase):
    value: str


class GetSubstanceByIdentifierQueryResult(QueryResultBase):
    substance: Substance


class GetSubstanceByIdentifierQueryHandler(
    QueryHandler[GetSubstanceByIdentifierQuery, GetSubstanceByIdentifierQueryResult]
):
    def __init__(self, unit_of_work: UnitOfWork) -> None:
        self._unit_of_work = unit_of_work

    def handle(
        self, query: GetSubstanceByIdentifierQuery
    ) -> GetSubstanceByIdentifierQueryResult:
        repo: SubstanceReadRepository = self._unit_of_work.get_repository(
            SubstanceReadRepository
        )  # type: ignore

        substance = repo.find_by_inchi_key(InchiKey(query.value))

        return GetSubstanceByIdentifierQueryResult(substance=substance)
