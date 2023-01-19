from shared.application.queries import QueryBase, QueryHandler, QueryResultBase
from shared.domain.uow import UnitOfWork
from substances_core.domain.syncs.enums import SyncState
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import SyncReadRepository


class GetSyncByIdQuery(QueryBase):
    """GetSyncByIdQuery"""

    sync_id: SyncId


class GetSyncByIdQueryResult(QueryResultBase):
    """GetSyncByIdQueryResult"""

    sync_id: SyncId
    provider: str
    state: SyncState


class GetSyncByIdQueryHandler(QueryHandler[GetSyncByIdQuery, GetSyncByIdQueryResult]):
    """GetSyncByIdQueryHandler"""

    def __init__(self, unit_of_work: UnitOfWork) -> None:
        self._unit_of_work = unit_of_work

    def handle(self, query: GetSyncByIdQuery) -> GetSyncByIdQueryResult:
        with self._unit_of_work as uow:
            repo: SyncReadRepository = uow.get_repository(SyncReadRepository)
            sync = repo.find_by_id(query.sync_id)
        return GetSyncByIdQueryResult(
            sync_id=sync.sync_id, provider=sync.provider, state=sync.state
        )
