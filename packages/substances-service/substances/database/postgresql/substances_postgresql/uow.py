import logging
import types
from typing import Any, Callable, Dict, List, Optional, Type

from shared.domain.events import DomainEvent, EventScope
from shared.domain.repositories import Repository
from shared.domain.uow import UnitOfWork
from shared.infrastructure.events import DomainEventDispatcher
from sqlalchemy.orm import Session
from substances_core.domain.susbtances.repository import SubstanceWriteRepository
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_postgresql.subtances.repository import (
    SubstancePostgresqlWriteRepository,
)
from substances_postgresql.syncs.repository import (
    SyncPostgresqlReadRepository,
    SyncPostgresqlWriteRepository,
)

LOGGER = logging.getLogger(__name__)


class PostgresqlUnitOfWork(UnitOfWork):
    """Unit of work implementation for postgresql with SQLAlchemy."""

    def __init__(
        self, session_factory: Callable[[], Any], dispatcher: DomainEventDispatcher
    ) -> None:
        self._session: Optional[Session] = None
        self._session_factory = session_factory
        self._dispatcher = dispatcher
        self._repositories: Dict[Type[Repository], Callable[[Session], Repository]] = {
            SyncReadRepository: SyncPostgresqlReadRepository,
            SyncWriteRepository: SyncPostgresqlWriteRepository,
            SubstanceWriteRepository: SubstancePostgresqlWriteRepository,
        }

    def __enter__(self):  # type: ignore
        if isinstance(self._session_factory, types.GeneratorType):
            self._session = next(self._session_factory)
        else:
            self._session = self._session_factory()
        return super().__enter__()

    @property
    def repositories(self) -> Dict[type, Callable[[Session], Repository]]:
        return self._repositories

    def get_repository(self, repository: type[Repository]) -> Repository:
        if repository not in self.repositories:
            raise Exception("No repository found!")

        if self._session is None:
            raise RuntimeError

        return self._repositories[repository](self._session)

    def commit(self, events: List[DomainEvent] = []) -> None:
        internal_events = filter(lambda ev: ev.scope is EventScope.INTERNAL, events)
        external_events = filter(lambda ev: ev.scope is EventScope.BOTH, events)

        # Execute domain side effects
        for ie in internal_events:
            self._dispatcher.dispatch(ie)

        # Save external notifications into message box

        # Persist changes
        if self._session:
            self._session.commit()

        # Dispatch external notification
        for ee in external_events:
            self._dispatcher.dispatch(ee)
