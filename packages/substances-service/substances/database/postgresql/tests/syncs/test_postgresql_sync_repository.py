import pytest
from shared.infrastructure.events import InternalDomainEventDispatcher
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_postgresql.syncs.entities import SyncItemModel, SyncModel
from substances_postgresql.syncs.repository import (
    SyncPostgresqlReadRepository,
    SyncPostgresqlWriteRepository,
)
from substances_postgresql.uow import PostgresqlUnitOfWork


@pytest.fixture()
def database_engine():
    engine = create_engine(
        "postgresql://localhost/substances_test?user=admin&password=admin"
    )
    return engine


@pytest.fixture()
def database_session():
    engine = create_engine(
        "postgresql://localhost/substances_test?user=admin&password=admin"
    )
    session = sessionmaker(autocommit=False, autoflush=False, bind=engine)
    yield session


def test_find_by_id(database_session):
    with PostgresqlUnitOfWork(
        database_session, dispatcher=InternalDomainEventDispatcher()
    ) as uow:
        sync_id = SyncId("4ab41039-f127-47fd-b073-482bb5086e73")
        read_repo: SyncPostgresqlReadRepository = uow.get_repository(SyncReadRepository)
        sync = read_repo.find_by_id(sync_id)

    assert sync is not None
