import logging
from typing import Optional
from uuid import UUID, uuid4

from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import Session
from sqlalchemy.orm.exc import NoResultFound
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.sync import Sync
from substances_postgresql.syncs.entities import SyncItemModel, SyncModel

LOGGER = logging.getLogger(__name__)


class SyncPostgresqlReadRepository(SyncReadRepository):
    """
    SyncPostgresqlReadRepository
    """

    def __init__(self, session: Session) -> None:
        self._session = session

    def find_by_id(self, sync_id: SyncId) -> Optional[Sync]:
        """
        Returns a 'Sync' aggregate root with the given identifier.
        """
        model: SyncModel = (
            self._session.query(SyncModel).filter_by(_sync_id=UUID(str(sync_id))).one()
        )

        return Sync(
            sync_id=model.sync_id,
            provider=model.provider,
            data_location=model.data_location,
            state=model.state,
            items=[
                SyncItem(
                    token=item.token,
                    location=item.location,
                    name=item.name,
                    transformation=item.transformation,
                    status=item.status,
                )
                for item in model.items
            ],
            created_at=model.created_at,
        )

    def find_last_one_by_provider(self, provider: str) -> Optional[Sync]:
        """
        Returns a reference to the last created sync for a specific provider.
        """
        try:
            model: SyncModel = (
                self._session.query(SyncModel).filter_by(provider=provider).one()
            )
            return Sync(
                sync_id=model.sync_id,
                provider=model.provider,
                state=model.state,
            )
        except NoResultFound as err:
            LOGGER.error(err)
        except Exception as err:
            LOGGER.error(err)
        finally:
            return None


class SyncPostgresqlWriteRepository(SyncWriteRepository):
    """
    SyncPostgresqlWriteRepository
    """

    def __init__(self, session: Session) -> None:
        self._session = session

    def next_id(self) -> SyncId:
        """next_id"""
        return SyncId(str(uuid4()))

    def add(self, sync: Sync) -> None:
        """
        Add a new 'Sync' into the persistent store.
        """
        sync_model = SyncModel(
            sync_id=sync.sync_id,
            provider=sync.provider,
            state=sync.state,
            data_location=sync.data_location,
        )
        self._session.add(sync_model)

    def update(self, sync: Sync) -> None:
        self._session.query(SyncModel).filter(SyncModel.sync_id == sync.sync_id).update(
            {SyncModel.state: sync.state, SyncModel.data_location: sync.data_location}
        )

    def insert_or_update_items(self, sync: Sync) -> None:
        for item in sync.items:
            self._session.execute(
                insert(SyncItemModel)
                .values(SyncItemModel.dict_from_domain(sync.sync_id, item))
                .on_conflict_do_update(
                    constraint=SyncItemModel.__table__.primary_key,
                    set_={
                        "name": item.name,
                        "transformation": item.transformation,
                        "status": item.status,
                    },
                )
            )
        self.update(sync)
