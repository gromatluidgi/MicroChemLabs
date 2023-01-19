import inspect
import logging
from typing import Callable, Dict, List, Optional, Type

from shared.application.events import DomainEventDispatcher
from shared.domain.events import DomainEvent
from shared.domain.repositories import Repository
from shared.domain.uow import UnitOfWork
from substances_core.domain.syncs.repository import SyncReadRepository
from substances_inmemory.context import MemoryContext
from substances_inmemory.syncs.repository import SyncInMemoryReadRepository

LOGGER = logging.getLogger(__name__)


class MemoryUnitOfWork(UnitOfWork):
    """Memory unit of work used for testing purpose."""

    def __init__(
        self,
        context: MemoryContext,
        dispatcher: Optional[DomainEventDispatcher[DomainEvent]] = None,
    ) -> None:
        self._context = context
        self._dispatcher = dispatcher
        self._repositories: Dict[
            Type[Repository], Callable[[MemoryContext], Repository]
        ] = {SyncReadRepository: SyncInMemoryReadRepository}

    @property
    def repositories(self) -> Dict[type, Callable[[MemoryContext], Repository]]:
        return self._repositories

    def get_repository(self, repository: type[Repository]) -> Repository:
        if not inspect.isabstract(repository):
            raise TypeError(
                "An abstract Repository type must be provided for retrieve a repository."
            )

        if repository not in self.repositories:
            raise Exception("No repository found!")
        return self._repositories[repository](self._context)

    def commit(self, events: List[DomainEvent] = []) -> None:
        LOGGER.debug("No transaction management for memory context.")
        if self._dispatcher is not None:
            pass
