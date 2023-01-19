from dataclasses import dataclass
from uuid import UUID

SyncItemId = UUID


@dataclass(frozen=True, eq=True)
class SyncId:

    identity: str

    def __str__(self) -> str:
        return self.identity
