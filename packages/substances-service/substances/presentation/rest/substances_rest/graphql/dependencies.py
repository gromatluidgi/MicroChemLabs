from functools import lru_cache

from shared.infrastructure.events import InternalDomainEventDispatcher
from sqlalchemy import create_engine
from sqlalchemy.orm import Session, scoped_session, sessionmaker
from substances_core.application.syncs.queries.get_sync_by_id import (
    GetSyncByIdQueryHandler,
)
from substances_postgresql.uow import PostgresqlUnitOfWork

from ..config import settings


# Database
@lru_cache()
def get_database_engine():
    return create_engine(str(settings.DB_URI))


def postgres_context(engine=get_database_engine):
    session: Session = scoped_session(
        sessionmaker(autocommit=False, autoflush=False, bind=engine())
    )
    try:
        yield session
    finally:
        session.close()


def domain_event_dispatcher() -> InternalDomainEventDispatcher:
    dispatcher = InternalDomainEventDispatcher()
    return dispatcher


# Transaction Management
def postgres_unit_of_work(
    session_factory=postgres_context, dispatcher=domain_event_dispatcher
):

    return PostgresqlUnitOfWork(session_factory(), dispatcher())


# Queries
def get_sync_by_id_query_handler(
    uow=postgres_unit_of_work,
) -> GetSyncByIdQueryHandler:
    return GetSyncByIdQueryHandler(uow())
