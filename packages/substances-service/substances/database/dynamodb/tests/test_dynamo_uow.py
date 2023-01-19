from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_dynamodb.uow import DynamoUnitOfWork


def test_uow():
    uow = DynamoUnitOfWork({})
    with uow as unit_of_work:
        print(unit_of_work)


def test_get_repository():
    table_names = {
        SyncReadRepository: "",
        SyncWriteRepository: "",
        SubstanceReadRepository: "",
        SubstanceWriteRepository: "",
    }
    uow = DynamoUnitOfWork(table_names)

    sync_read_repository: SyncReadRepository = uow.get_repository(SyncReadRepository)

    assert sync_read_repository is not None
