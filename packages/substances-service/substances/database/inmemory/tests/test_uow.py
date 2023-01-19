# from core.domain.syncs.repository import SyncReadRepository
# from database.inmemory.context import MemoryContext
# from database.inmemory.syncs.repository import SyncInMemoryReadRepository
# from database.inmemory.uow import MemoryUnitOfWork


# def test_get_repository():
#     context = MemoryContext()
#     uow = MemoryUnitOfWork(context)

#     # Act
#     repository = uow.get_repository(SyncReadRepository)

#     # Assert
#     assert repository is not None
#     assert isinstance(repository, SyncInMemoryReadRepository)
