from abc import ABC, abstractmethod
from typing import Generic, List

from shared.domain.aggregates import TAggregateRoot
from substances_core.domain.syncs.entities.item import SyncItem


class SyncLoader(ABC, Generic[TAggregateRoot]):
    """
    SyncLoader
    """

    @abstractmethod
    def load(self, item: SyncItem) -> List[TAggregateRoot]:
        raise NotImplementedError


class SyncLoaderFactory(ABC):
    """SyncLoaderFactory"""

    @abstractmethod
    def from_provider(self, provider: str) -> SyncLoader[TAggregateRoot]:
        raise NotImplementedError
