import logging
import uuid

from shared.application.events import DomainEventHandler
from shared.domain.uow import UnitOfWork
from substances_core.domain.susbtances.repository import SubstanceWriteRepository
from substances_core.domain.syncs.enums import SyncItemStatus, SyncState
from substances_core.domain.syncs.events import SyncCheckedOutEvent
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.services.loader import SyncLoaderFactory
from substances_core.domain.syncs.sync import Sync

LOGGER = logging.getLogger()

# WARNING: only for demo purpose


class SyncCheckoutEventHandler(DomainEventHandler[SyncCheckedOutEvent]):
    """SyncCheckedOutEventHandler"""

    def __init__(
        self, unit_of_work: UnitOfWork, sync_loader_factory: SyncLoaderFactory
    ) -> None:
        self._unit_of_work = unit_of_work
        self._sync_loader_factory = sync_loader_factory

    def handle(self, event: SyncCheckedOutEvent) -> None:
        LOGGER.debug("SyncCheckedOutEvent Handled")

        # . Fetch sync
        sync_read_repo: SyncReadRepository = self._unit_of_work.get_repository(
            SyncReadRepository
        )
        sync = sync_read_repo.find_by_id(uuid.UUID(event.aggregate_id))

        if sync is None:
            raise ValueError("Sync with id not found.")

        # . State validation
        self.__validate_state(sync)

        # . Execute business logic
        loader = self._sync_loader_factory.from_provider(sync.provider)

        substance_write_repo: SubstanceWriteRepository = (
            self._unit_of_work.get_repository(SubstanceWriteRepository)
        )

        # . Load transformation from each SyncItem
        for item in sync.items:
            substance_write_repo.batch_insert_or_update(loader.load(item))
            item.status = SyncItemStatus.LOADED

        try:
            sync.complete()
        except RuntimeError:
            sync.state = SyncState.ERROR

        # Update sync
        sync_write_repo: SyncWriteRepository = self._unit_of_work.get_repository(
            SyncWriteRepository
        )
        sync_write_repo.insert_or_update_items(sync)

        # . End of transaction: persist changes
        self._unit_of_work.commit()

    def __validate_state(self, sync: Sync):
        if sync.state is not SyncState.CHECKOUT:
            raise RuntimeError()

        if len(sync.items) == 0:
            raise RuntimeError("Sync has no SyncItem to load into database.")
