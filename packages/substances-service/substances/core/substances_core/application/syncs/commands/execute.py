import logging
from typing import Optional, Set

from shared.application.commands import CommandBase, CommandHandler, CommandResultBase
from shared.domain.uow import UnitOfWork
from shared.infrastructure.queue import Message, MessageQueue, MessageType
from substances_core.application.syncs.commands.internal.extract import (
    ExtractSyncCommand,
)
from substances_core.domain.syncs.enums import SyncState
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.sync import Sync

LOGGER = logging.getLogger(__name__)


class ExecuteSyncCommand(CommandBase):
    """ExecuteSyncCommand"""

    provider: str
    """Indicates where to fetch the data."""

    indexes: Set[str] = set()
    """Indexes can be provided to extract specific datasets."""

    data_location: Optional[str] = None


class ExecuteSyncCommandResult(CommandResultBase):
    """ScheduleSyncResult"""

    sync_id: str
    """Unique id used to track the state of a synchronization process."""

    state: SyncState
    """Current state of the sync machine."""


class ExecuteSyncCommandHandler(
    CommandHandler[ExecuteSyncCommand, ExecuteSyncCommandResult]
):
    """ExecuteSyncCommandHandler"""

    def __init__(
        self,
        unit_of_work: UnitOfWork,
        message_queue: Optional[MessageQueue] = None,
    ) -> None:
        self._unit_of_work = unit_of_work
        self._message_queue = message_queue

    def execute(self, command: ExecuteSyncCommand) -> ExecuteSyncCommandResult:
        """Execute use case."""

        LOGGER.info("Prepare ETL execution process.")

        with self._unit_of_work as uow:

            # . Retrieve repository from uow
            sync_read_repository: SyncReadRepository = uow.get_repository(
                SyncReadRepository
            )  # type: ignore
            sync_write_repository: SyncWriteRepository = uow.get_repository(
                SyncWriteRepository
            )  # type: ignore

            # . Fetch last sync for provider
            last_sync = sync_read_repository.find_last_one_by_provider(command.provider)

            # . Business rules checking
            if last_sync and not self._can_schedule_new_sync(last_sync):
                raise RuntimeError(f"Can't schedule new sync for ${command.provider}")

            # . Create new sync
            sync_id = sync_write_repository.next_id()
            sync = Sync(
                sync_id=sync_id,
                provider=command.provider,
                data_location=command.data_location,
                indexes=command.indexes,
            )
            sync_write_repository.add(sync)

            # . Commit transaction
            uow.commit(events=sync.domain_events)

        # . (Background Thread) Send Sync aggregate to queue worker for futher processing
        # PoF: compensated by sync state business validation + outbox worker
        try:
            if self._message_queue is not None:
                self._message_queue.add_message(
                    Message(
                        type=MessageType.COMMAND,
                        body=ExtractSyncCommand(
                            sync_id=str(sync_id), indexes=command.indexes
                        ),
                    )
                )
                self._message_queue.execute()
            else:
                LOGGER.debug("No message queue provided for ExecuteSync command.")
        except Exception as err:  # pylint: disable=broad-except
            LOGGER.error(err)

        return ExecuteSyncCommandResult(sync_id=str(sync_id), state=sync.state)

    @classmethod
    def _can_schedule_new_sync(cls, last_sync: Sync) -> bool:
        return last_sync.is_synced()
