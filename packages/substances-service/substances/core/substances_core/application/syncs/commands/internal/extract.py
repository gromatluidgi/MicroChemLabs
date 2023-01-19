import logging
from typing import Optional

from shared.application.commands import CommandBase, CommandHandler, CommandResultBase
from shared.domain.uow import UnitOfWork
from shared.infrastructure.queue import Message, MessageQueue, MessageType
from substances_core.application.syncs.commands.internal.transform import (
    TransformSyncCommand,
)
from substances_core.domain.syncs.enums import SyncState
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.services.extractor import SyncExtractorFactory
from substances_core.domain.syncs.sync import Sync

LOGGER = logging.getLogger(__name__)


class ExtractSyncCommand(CommandBase):
    """ExecuteSyncCommand"""

    sync_id: str


class ExtractSyncCommandResult(CommandResultBase):
    """ScheduleSyncResult"""

    sync_id: str
    state: SyncState
    has_items: bool


class ExtractSyncCommandHandler(
    CommandHandler[ExtractSyncCommand, ExtractSyncCommandResult]
):
    """ExtractSyncCommandHandler"""

    def __init__(
        self,
        unit_of_work: UnitOfWork,
        sync_extractor_factory: SyncExtractorFactory,
        message_queue: Optional[MessageQueue] = None,
    ) -> None:
        self._unit_of_work = unit_of_work
        self._sync_extractor_factory = sync_extractor_factory
        self._message_queue = message_queue

    def execute(self, command: ExtractSyncCommand) -> ExtractSyncCommandResult:
        """Execute use case."""

        LOGGER.info("ETL - Start extraction process")

        with self._unit_of_work as uow:
            # . Retrieve repository from uow
            sync_read_repository: SyncReadRepository = uow.get_repository(
                SyncReadRepository
            )  # type: ignore
            sync_write_repository: SyncWriteRepository = uow.get_repository(
                SyncWriteRepository
            )  # type: ignore

            # . Fetch sync
            sync = sync_read_repository.find_by_id(SyncId(command.sync_id))

            if sync is None:
                raise ValueError("Sync with id not found.")

            # Execute business logic
            extractor = self._sync_extractor_factory.from_provider(
                provider=sync.provider,
                data_location=sync.data_location,
                indexes=sync.indexes,
            )
            sync_items = extractor.extract()

            # . Update sync and sync items
            try:
                sync.checkin(
                    data_location=extractor.get_data_location(), items=sync_items
                )
                sync_write_repository.insert_or_update_items(sync)
            except Exception as err:
                LOGGER.error(err)
                raise

            # . Commit transaction
            uow.commit(sync.domain_events)

        # . Send an async checkout command
        # Eventual Consistency:
        try:
            if self._message_queue is not None:
                self._message_queue.add_message(
                    Message(
                        type=MessageType.COMMAND,
                        body=TransformSyncCommand(sync_id=command.sync_id),
                    )
                )
                self._message_queue.execute()
        except Exception as err:
            LOGGER.error(err)

        return ExtractSyncCommandResult(
            sync_id=str(sync.sync_id), state=sync.state, has_items=sync.has_items()
        )

    def _can_schedule_new_sync(self, last_sync: Sync) -> bool:
        return last_sync.is_synced()
