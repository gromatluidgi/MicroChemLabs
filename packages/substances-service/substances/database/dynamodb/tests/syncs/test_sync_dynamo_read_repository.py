import datetime

from moto import mock_dynamodb
from substances_core.domain.syncs.sync import Sync
from substances_dynamodb.context import DynamoContext
from substances_dynamodb.syncs.repository import (
    SyncDynamoReadRepository,
    SyncDynamoWriteRepository,
)


@mock_dynamodb
def test_find_by_id():
    # Prepare
    context = DynamoContext()
    context._create_syncs_table()
    read_repository = SyncDynamoReadRepository("syncs")
    write_repository = SyncDynamoWriteRepository("syncs")
    sync = Sync(sync_id=write_repository.next_id(), provider="PUBCHEM")
    write_repository.add(sync)

    # Act
    item = read_repository.find_by_id(sync.sync_id)

    # Assert
    assert item is not None
    assert item.sync_id == sync.sync_id
    assert item.created_at == sync.created_at.replace(microsecond=0)


@mock_dynamodb
def test_find_last_one_by_provider():
    # Prepare
    context = DynamoContext()
    context._create_syncs_table()
    read_repository = SyncDynamoReadRepository("syncs")
    write_repository = SyncDynamoWriteRepository("syncs")
    sync_one = Sync(write_repository.next_id(), "PUBCHEM", data_location="test_one")
    sync_two = Sync(write_repository.next_id(), "PUBCHEM", data_location="test_two")
    sync_three = Sync(
        sync_id=write_repository.next_id(),
        provider="PUBCHEM",
        data_location="test_three",
        created_at=datetime.datetime.now().replace(day=1, minute=1, second=1),
    )
    write_repository.add(sync_one)
    write_repository.add(sync_two)
    write_repository.add(sync_three)

    # Act
    item = read_repository.find_last_one_by_provider("PUBCHEM")

    # Assert
    assert item is not None
    assert item.sync_id == sync_two.sync_id
    assert item.created_at == sync_two.created_at.replace(microsecond=0)
