# from core.domain.syncs.enums import SyncState
# from core.domain.syncs.sync import Sync
# from database.inmemory.context import MemoryContext
# from database.inmemory.syncs.repository import SyncInMemoryWriteRepository


# def test_add():
#     context = MemoryContext()
#     repository = SyncInMemoryWriteRepository(context)

#     repository.add(Sync("test", "test", SyncState.PENDING))

#     assert len(context.syncs) > 0
