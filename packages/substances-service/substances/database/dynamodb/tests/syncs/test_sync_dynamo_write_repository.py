import datetime

from moto import mock_dynamodb
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus, SyncState
from substances_core.domain.syncs.sync import Sync
from substances_dynamodb.context import DynamoContext
from substances_dynamodb.syncs.repository import (
    SyncDynamoReadRepository,
    SyncDynamoWriteRepository,
)


@mock_dynamodb
def test_add():
    # Prepare
    context = DynamoContext()
    context._create_syncs_table()
    write_repository = SyncDynamoWriteRepository("syncs")
    sync = Sync(
        sync_id=write_repository.next_id(),
        provider="PUBCHEM",
    )

    # Act
    write_repository.add(sync)

    # Assert
    table = write_repository.get_table()
    assert table.item_count == 1


@mock_dynamodb
def test_update():
    # Prepare
    context = DynamoContext()
    context._create_syncs_table()
    read_repository = SyncDynamoReadRepository("syncs")
    write_repository = SyncDynamoWriteRepository("syncs")
    sync = Sync(
        id=write_repository.next_id(),
        provider="PUBCHEM",
        updated_at=datetime.datetime.now(),
    )
    write_repository.add(sync)

    # Act
    sync.state = SyncState.COMPLETED
    write_repository.update(sync)
    item = read_repository.find_by_id(sync.sync_id)

    # Assert
    assert item is not None
    assert item.state == SyncState.COMPLETED


@mock_dynamodb
def test_insert_or_update_items():
    # Prepare
    context = DynamoContext()
    context._create_syncs_table()
    read_repository = SyncDynamoReadRepository("syncs")
    write_repository = SyncDynamoWriteRepository("syncs")
    sync = Sync(
        sync_id=write_repository.next_id(),
        provider="PUBCHEM",
    )
    write_repository.add(sync)

    # Act
    sync.items = [
        SyncItem(token="token1", location="location1"),
        SyncItem(token="token2", location="location2"),
    ]
    write_repository.insert_or_update_items(sync)

    sync.items = [
        SyncItem(
            token="token1",
            location="location1",
            status=SyncItemStatus.ERROR,
        ),
        SyncItem(token="token2", location="location2"),
    ]
    write_repository.insert_or_update_items(sync)

    item = read_repository.find_by_id(sync.sync_id)

    # Assert
    assert item is not None
    assert len(item.items) > 0
