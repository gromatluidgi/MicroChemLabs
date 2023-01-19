import uuid
from typing import Optional

from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.sync import Sync
from substances_inmemory.context import MemoryContext


class SyncInMemoryReadRepository(SyncReadRepository):
    """
    SyncInMemoryReadRepository
    """

    def __init__(self, context: MemoryContext) -> None:
        self._syncs = context.syncs
        self._context = context

    def find_by_id(self, sync_id: SyncId) -> Optional[Sync]:
        """
        Returns a reference to an 'Sync' aggregate root with the given identifier.
        """
        return self._syncs.get(sync_id)

    def find_last_one_by_provider(self, provider: str) -> Optional[Sync]:
        """
        Returns a reference to the last created sync for a specific provider.
        """
        filtered_syncs = filter(
            lambda s: s.provider.casefold() == provider.casefold(), self._syncs.values()
        )
        ordered_syncs = sorted(filtered_syncs, key=lambda s: s.created_at, reverse=True)
        try:
            return ordered_syncs.pop(0)
        except IndexError:
            return None


class SyncInMemoryWriteRepository(SyncWriteRepository):
    """
    SyncInMemoryWriteRepository
    """

    def __init__(self, context: MemoryContext) -> None:
        self._syncs = context.syncs
        self._context = context

    def next_id(self) -> str:
        """next_id"""
        return str(uuid.uuid4())

    def add(self, sync: Sync) -> None:
        """
        Add a new 'Sync' into the persistent store.
        """
        if sync.sync_id in list(self._syncs):
            raise KeyError(f"Sync {sync.sync_id} already in the collection.")
        self._syncs[sync.sync_id] = sync

    def insert_or_update_items(self, sync: Sync) -> None:
        return super().insert_or_update_items(sync)
