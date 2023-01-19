import logging

from shared.application.commands import CommandBase, CommandHandler, CommandResultBase
from shared.domain.uow import UnitOfWork
from substances_core.domain.susbtances.repository import SubstanceWriteRepository
from substances_core.domain.susbtances.substance import Substance
from substances_core.domain.syncs.enums import SyncState
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.services.loader import SyncLoader, SyncLoaderFactory
from substances_core.domain.syncs.sync import Sync

LOGGER = logging.getLogger()


class LoadSyncCommand(CommandBase):
    """LoadSyncCommand"""

    sync_id: str


class LoadSyncCommandResult(CommandResultBase):
    """LoadSyncCommandResult"""

    sync_id: str
    state: SyncState
    items_loaded: int


class LoadSyncCommandHandler(CommandHandler[LoadSyncCommand, LoadSyncCommandResult]):
    """LoadSyncCommandHandler"""

    def __init__(
        self, unit_of_work: UnitOfWork, sync_loader_factory: SyncLoaderFactory
    ) -> None:
        self._unit_of_work = unit_of_work
        self._sync_loader_factory = sync_loader_factory

    def execute(self, command: LoadSyncCommand) -> LoadSyncCommandResult:
        LOGGER.info("ETL - Start loading process")

        # 1. Open transaction
        with self._unit_of_work as uow:
            # 2. Fetch sync
            sync_read_repo: SyncReadRepository = uow.get_repository(SyncReadRepository)  # type: ignore
            sync = sync_read_repo.find_by_id(SyncId(command.sync_id))

            if sync is None:
                raise ValueError("Sync with id not found.")

            # 3. Current state validation
            self.__validate_state(sync)

            # . Load transformation from each SyncItem
            loader: SyncLoader[Substance] = self._sync_loader_factory.from_provider(
                sync.provider
            )
            substance_write_repo: SubstanceWriteRepository = uow.get_repository(
                SubstanceWriteRepository
            )  # type: ignore

            substances = []
            items_loaded = 0
            for item in sync.items:
                try:
                    substances = [*substances, *loader.load(item)]
                    sync.item_loaded(item.token)
                except Exception as err:
                    LOGGER.error(err)
                    sync.item_load_failed(item.token)

            # TODO: can be improved
            try:
                sync.complete()  # if a SyncItem isn't loaded, an exception is raised
                substance_write_repo.batch_insert_or_update(
                    substances
                )  # eventual infrastructure exception
            except Exception as err:
                LOGGER.error(err)
                sync.failed()
            else:
                items_loaded = len(
                    substances
                )  # TODO: not accurate -> batch_insert can have unprocessed items

            # Update sync
            sync_write_repo: SyncWriteRepository = uow.get_repository(
                SyncWriteRepository
            )  # type: ignore
            sync_write_repo.insert_or_update_items(sync)

            # 5. End of transaction
            uow.commit()

        return LoadSyncCommandResult(
            sync_id=str(sync.sync_id), state=sync.state, items_loaded=items_loaded
        )

    def __validate_state(self, sync: Sync) -> None:
        if sync.state is not SyncState.CHECKOUT:
            raise RuntimeError()
