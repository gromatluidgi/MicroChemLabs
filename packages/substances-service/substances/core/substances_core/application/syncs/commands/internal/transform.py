import logging
from typing import Optional

from shared.application.commands import CommandBase, CommandHandler, CommandResultBase
from shared.domain.aggregates import AggregateRoot
from shared.domain.uow import UnitOfWork
from shared.infrastructure.queue import Message, MessageQueue, MessageType
from substances_core.application.syncs.commands.internal.load import LoadSyncCommand
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.services.transformer import (
    SyncTransformer,
    SyncTransformerFactory,
)

LOGGER = logging.getLogger(__name__)


class TransformSyncCommand(CommandBase):
    """TransformSyncCommand"""

    sync_id: str


class TransformSyncCommandResult(CommandResultBase):
    """TransformSyncCommandResult"""

    sync_id: str


class TransformSyncCommandHandler(
    CommandHandler[TransformSyncCommand, TransformSyncCommandResult]
):
    """TransformSyncCommandHandler"""

    def __init__(
        self,
        unit_of_work: UnitOfWork,
        sync_transformer_factory: SyncTransformerFactory,
        message_queue: Optional[MessageQueue] = None,
    ) -> None:
        self._unit_of_work = unit_of_work
        self._sync_transformer_factory = sync_transformer_factory
        self._message_queue = message_queue

    def execute(self, command: TransformSyncCommand) -> TransformSyncCommandResult:
        """Execute use case."""

        LOGGER.info("ETL - Start transformation process")

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

            # . Validate and transform sync item payload into domain model (binary form)
            transformer: SyncTransformer[
                AggregateRoot
            ] = self._sync_transformer_factory.from_provider(
                sync.provider
            )  # TODO: review

            for item in sync.items:
                try:
                    transformation = transformer.transform(item)
                    sync.item_transformed(item.token, transformation)
                except Exception as err:
                    LOGGER.error(err)
                    sync.item_transformation_failed(item.token)

            # . Update sync and sync items
            try:
                sync.checkout()
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
                        body=LoadSyncCommand(sync_id=command.sync_id),
                    )
                )
                self._message_queue.execute()
        except Exception as err:
            LOGGER.error(err)

        return TransformSyncCommandResult(sync_id=str(sync.sync_id))
