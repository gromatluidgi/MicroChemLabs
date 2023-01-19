from typing import Dict

from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.sync import Sync


class MemoryContext:
    def __init__(self) -> None:
        self._syncs: Dict[SyncId, Sync] = {}

    @property
    def syncs(self) -> Dict[SyncId, Sync]:
        return self._syncs
