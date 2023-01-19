import pytest
from shared.infrastructure.events import InternalDomainEventDispatcher
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus, SyncState
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.sync import Sync
from substances_postgresql.subtances.entities import SubstanceModel
from substances_postgresql.syncs.entities import SyncItemModel, SyncModel
from substances_postgresql.syncs.repository import (
    SyncPostgresqlReadRepository,
    SyncPostgresqlWriteRepository,
)
from substances_postgresql.uow import PostgresqlUnitOfWork


@pytest.fixture()
def database_engine():
    engine = create_engine("postgresql://admin:admin@localhost/substances_test")
    return engine


@pytest.fixture()
def database_session():
    engine = create_engine(
        "postgresql://localhost/substances_test?user=admin&password=admin"
    )
    session = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    try:
        yield session
    finally:
        session.close()


def test_remove_all_tables(database_engine):
    SyncItemModel.__table__.drop(database_engine)
    SyncModel.__table__.drop(database_engine)
    SubstanceModel.__table__.drop(database_engine)


def test_create_all_tables(database_engine):
    SyncModel.__table__.create(bind=database_engine, checkfirst=True)
    SyncItemModel.__table__.create(bind=database_engine, checkfirst=True)
    SubstanceModel.__table__.create(bind=database_engine, checkfirst=True)


def test_commit(database_session):
    # Arrange
    with PostgresqlUnitOfWork(
        database_session, dispatcher=InternalDomainEventDispatcher()
    ) as uow:
        read_repo: SyncPostgresqlReadRepository = uow.get_repository(SyncReadRepository)
        write_repo: SyncPostgresqlWriteRepository = uow.get_repository(
            SyncWriteRepository
        )

        sync_id = write_repo.next_id()
        sync = Sync(sync_id, "test", SyncState.PENDING)
        sync.items.append(SyncItem("test", "test", SyncItemStatus.PENDING))
        write_repo.add(sync)

        # Act
        uow.commit()
        sync = read_repo.find_by_id(sync_id)

    # Assert
    assert sync is not None
