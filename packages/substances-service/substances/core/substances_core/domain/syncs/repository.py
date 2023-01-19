from abc import ABC, abstractmethod
from typing import Optional

from shared.domain.repositories import Repository
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.sync import Sync


class SyncReadRepository(Repository, ABC):
    """SyncReadRepository"""

    @abstractmethod
    def find_by_id(self, sync_id: SyncId) -> Optional[Sync]:
        """
        Returns a reference to an 'Sync' aggregate root with the given identifier.
        """
        raise NotImplementedError()

    @abstractmethod
    def find_last_one_by_provider(self, provider: str) -> Optional[Sync]:
        """
        Returns a reference to the last created sync for a specific provider.
        """
        raise NotImplementedError()


class SyncWriteRepository(Repository, ABC):
    """SyncWriteRepository"""

    @abstractmethod
    def next_id(self) -> SyncId:
        """next_id"""
        raise NotImplementedError()

    @abstractmethod
    def add(self, sync: Sync) -> None:
        """
        Add a new 'Sync' into the persistent store.
        """
        raise NotImplementedError()

    @abstractmethod
    def update(self, sync: Sync) -> None:
        """
        Update sync.
        """
        raise NotImplementedError()

    @abstractmethod
    def insert_or_update_items(self, sync: Sync) -> None:
        """
        Insert or update items for a Sync aggregate.
        """
        raise NotImplementedError()
