from abc import ABC, abstractmethod
from typing import Any, Generic, List, Union

from shared.domain.aggregates import TAggregateRoot
from substances_core.domain.syncs.entities.item import SyncItem


class SyncTransformer(ABC, Generic[TAggregateRoot]):
    """
    SyncTransforming
    """

    @abstractmethod
    def transform(self, item: SyncItem, inmemory: bool = False) -> Union[str, bytes]:
        raise NotImplementedError

    @abstractmethod
    def parse(self, data: Any, chunk_start=0, chunk_size=100) -> List[TAggregateRoot]:
        raise NotImplementedError


class SyncTransformerFactory(ABC):
    """Abstract domain service factory."""

    @abstractmethod
    def from_provider(self, provider: str) -> SyncTransformer[TAggregateRoot]:
        raise NotImplementedError
