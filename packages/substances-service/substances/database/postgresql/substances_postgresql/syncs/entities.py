"""SyncModel"""
import uuid
from typing import Any, Dict

from sqlalchemy import Column, Enum, ForeignKey, LargeBinary, String
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import declarative_base, relationship
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus, SyncState
from substances_core.domain.syncs.objects import SyncId
from substances_postgresql.mixins import TimestampMixin

Base = declarative_base()


class SyncModel(Base, TimestampMixin):  # type: ignore
    """SyncModel"""

    __tablename__ = "sync"

    _sync_id = Column("id", UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)
    provider = Column(String, index=True, nullable=False)
    data_location = Column(String, index=False, nullable=True)
    state = Column(
        Enum(SyncState, values_callable=lambda obj: [e.value for e in obj]),
        index=True,
        nullable=False,
    )
    items = relationship("SyncItemModel", backref="sync", lazy="immediate")  # type: ignore

    @hybrid_property
    def sync_id(self) -> SyncId:
        return SyncId(str(self._sync_id))

    @sync_id.setter
    def sync_id(self, value: SyncId) -> None:
        self._sync_id = value.identity


class SyncItemModel(Base, TimestampMixin):  # type: ignore
    """Table for storing sync items."""

    __tablename__ = "sync_item"
    token = Column(String, primary_key=True)
    location = Column(String, index=True, nullable=False)
    status = Column(
        Enum(SyncItemStatus, values_callable=lambda obj: [e.value for e in obj]),
        index=True,
        nullable=False,
    )
    name = Column(String, index=False, nullable=True)
    transformation = Column(LargeBinary, index=False, nullable=True)
    sync_id = Column(UUID(as_uuid=True), ForeignKey("sync.id"))

    @classmethod
    def from_domain(cls, sync_id: SyncId, item: SyncItem) -> "SyncItemModel":
        entity = cls()
        entity.sync_id = str(sync_id)
        entity.location = item.location
        entity.token = item.token
        entity.name = item.name
        entity.transformation = item.transformation
        entity.status = item.status
        return entity

    @classmethod
    def dict_from_domain(cls, sync_id: SyncId, item: SyncItem) -> Dict[str, Any]:
        return {
            "token": item.token,
            "location": item.location,
            "name": item.name,
            "transformation": item.transformation,
            "status": item.status,
            "sync_id": str(sync_id),
        }
