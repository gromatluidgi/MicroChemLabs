from abc import ABC, abstractmethod
from typing import Any, Callable, Dict, List, Type

from shared.domain.events import DomainEvent
from shared.domain.repositories import Repository


class UnitOfWork(ABC):
    """
    Unit of work pattern wrap database calls
    into a single transaction unit.
    """

    def __enter__(self) -> "UnitOfWork":
        return self

    def __exit__(self, *args) -> None:  # type: ignore
        pass

    @property
    @abstractmethod
    def repositories(self) -> Dict[Type[Repository], Callable[[Any], Repository]]:
        """Get a catalog of available repositories."""
        raise NotImplementedError

    @abstractmethod
    def get_repository(self, repository: Type[Repository]) -> Repository:
        """Retrieve a repository for an entity."""
        raise NotImplementedError

    @abstractmethod
    def commit(self, events: List[DomainEvent] = []) -> None:
        """Commit any changes occured during transaction."""
        raise NotImplementedError
