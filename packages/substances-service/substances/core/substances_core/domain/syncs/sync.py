import datetime
from dataclasses import dataclass, field
from typing import List, Optional, Set

from shared.domain.aggregates import AggregateRoot
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus, SyncState
from substances_core.domain.syncs.events import SyncCheckedInEvent, SyncExecutedEvent
from substances_core.domain.syncs.objects import SyncId


@dataclass
class Sync(AggregateRoot):
    """
    Synchronization
    """

    sync_id: SyncId
    provider: str
    state: SyncState = SyncState.PENDING
    indexes: Set[str] = field(default_factory=set)
    items: List[SyncItem] = field(default_factory=list)
    """List of items to process by ETL."""
    data_location: Optional[str] = field(default=None)
    """Path to the data that's need to be fetched from provider."""
    created_at: Optional[datetime.datetime] = field(default=None)
    """Creation date."""
    updated_at: Optional[datetime.datetime] = field(default=None)

    def __post_init__(self) -> None:
        """Post initialization logic."""
        if self.created_at is None:
            self.created_at = datetime.datetime.now(datetime.timezone.utc)
            self.add_domain_event(SyncExecutedEvent(str(self.sync_id), self.provider))

    def __can_checkin(self) -> None:
        if len(self.items) > 0:
            raise RuntimeError

    def checkin(self, data_location: str, items: List[SyncItem]) -> None:
        # Validation
        self.__can_checkin()

        # Mutation
        self.data_location = data_location

        if len(items) == 0:
            self.state = SyncState.COMPLETED
        else:
            for item in items:
                self.items.append(item)
            self.state = SyncState.CHECKIN
            # Add domain events
            self.add_domain_event(SyncCheckedInEvent(str(self.sync_id)))

        self.updated_at = datetime.datetime.now(datetime.timezone.utc)

    def __can_checkout(self) -> None:
        if len(self.items) == 0:
            return

        if (
            len(
                [
                    item
                    for item in self.items
                    if item.status != SyncItemStatus.TRANSFORMED
                ]
            )
            > 0
        ):
            raise RuntimeError("Can't checkout - An item has an invalid status.")

    def checkout(self) -> None:
        # Validation
        self.__can_checkout()

        # Mutation
        self.state = SyncState.CHECKOUT
        self.updated_at = datetime.datetime.now(datetime.timezone.utc)

        # Add domain events
        # self.add_domain_event(SyncCheckedOutEvent(str(self.id), SyncScope.SUBSTANCE))

    def __can_complete(self) -> bool:

        if (
            len([item for item in self.items if item.status != SyncItemStatus.LOADED])
            > 0
        ):
            return False
        return True

    def complete(self) -> None:
        self.updated_at = datetime.datetime.now(datetime.timezone.utc)
        if self.__can_complete():
            self.state = SyncState.COMPLETED
        else:
            raise RuntimeError()

    def has_items(self) -> bool:
        return len(self.items) > 0

    def get_pending_items(self) -> List[SyncItem]:
        return [item for item in self.items if item.status == SyncItemStatus.PENDING]

    def has_sync_error(self) -> bool:
        return (
            len([item for item in self.items if item.status == SyncItemStatus.ERROR])
            > 0
        )

    def is_synced(self) -> bool:
        return self.state == SyncState.COMPLETED or (
            len([item for item in self.items if item.status != SyncItemStatus.LOADED])
            == 0
        )

    def failed(self) -> None:
        self.state = SyncState.ERROR

    def item_transformed(self, token: str, transformation: bytes) -> None:
        item = next(it for it in self.items if it.token == token)
        item.update_transformed_data(transformation)
        item.status = SyncItemStatus.TRANSFORMED

    def item_transformation_failed(self, token: str) -> None:
        item = next(it for it in self.items if it.token == token)
        item.status = SyncItemStatus.ERROR

    def item_loaded(self, token: str) -> None:
        item = next(it for it in self.items if it.token == token)
        item.status = SyncItemStatus.LOADED

    def item_load_failed(self, token: str) -> None:
        item = next(it for it in self.items if it.token == token)
        item.status = SyncItemStatus.ERROR
