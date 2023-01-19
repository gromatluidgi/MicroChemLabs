import logging
from datetime import datetime, timezone
from typing import Optional
from uuid import uuid4

import boto3
from boto3.dynamodb.conditions import Key
from botocore.exceptions import ClientError
from substances_core.domain.syncs.enums import SyncState
from substances_core.domain.syncs.objects import SyncId
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_core.domain.syncs.sync import Sync
from substances_dynamodb.syncs.translator import SyncItemTranslator

LOGGER = logging.getLogger(__name__)


class SyncDynamoReadRepository(SyncReadRepository):
    """
    SyncDynamoReadRepository
    """

    def __init__(self, table_name: str) -> None:
        self._dyn_resource = boto3.resource("dynamodb")
        self._table = self._dyn_resource.Table(table_name)
        LOGGER.debug("Initialize Read Repository for table: ", table_name)

    def find_by_id(self, sync_id: SyncId) -> Optional[Sync]:
        """
        Gets sync data from the data store for a specific sync.
        """
        try:
            response = self._table.get_item(Key={"id": str(sync_id)})
        except ClientError as err:
            LOGGER.error(
                "Couldn't get sync %s from table %s. Here's why: %s: %s",
                sync_id,
                self._table.name,
                err.response["Error"]["Code"],
                err.response["Error"]["Message"],
            )
            return None
        else:
            item = response["Item"]
            return Sync(
                sync_id=SyncId(item["id"]),
                provider=item["provider"],
                state=SyncState(item["state"]),
                data_location=item["dataLocation"]
                if "dataLocation" in item.keys()
                else None,
                items=SyncItemTranslator.from_records(item["items"])
                if "items" in item.keys()
                else [],
                indexes=item["indexes"] if "indexes" in item.keys() else [],
                created_at=datetime.fromtimestamp(
                    int(item["createdAt"]), tz=timezone.utc
                ),
            )

    def find_last_one_by_provider(self, provider: str) -> Optional[Sync]:
        """
        Returns a reference to the last created sync for a specific provider.
        """
        try:
            response = self._table.query(
                IndexName="SyncProviderIndex",
                Limit=1,
                KeyConditionExpression=Key("provider").eq(provider),
                ScanIndexForward=False,
            )
        except Exception as err:
            LOGGER.error(err)
            return None
        else:
            items = response["Items"]
            if len(items) == 0:
                return None

            return Sync(
                sync_id=SyncId(items[0]["id"]),
                provider=items[0]["provider"],
                state=items[0]["state"],
                created_at=datetime.fromtimestamp(
                    items[0]["createdAt"], tz=timezone.utc
                ),
                data_location=items[0]["dataLocation"],
            )


class SyncDynamoWriteRepository(SyncWriteRepository):
    """
    SyncDynamoWriteRepository
    """

    def __init__(self, table_name: str) -> None:
        self._dyn_resource = boto3.resource("dynamodb")
        self._table = self._dyn_resource.Table(table_name)
        LOGGER.debug("Initialize Write Repository for table: ", table_name)

    def get_table(self):  # type: ignore
        return self._table

    def next_id(self) -> SyncId:
        """next_id"""
        return SyncId(str(uuid4()))

    def add(self, sync: Sync) -> None:
        """
        Adds a sync to the DynamoDB table.
        """
        self._table.put_item(
            Item={
                "id": str(sync.sync_id),
                "provider": sync.provider,
                "state": sync.state,
                "indexes": [index for index in sync.indexes],
                "dataLocation": sync.data_location,
                "createdAt": int(sync.created_at.timestamp()),
            }
        )

    def update(self, sync: Sync) -> None:
        self._table.update_item(
            Key={"id": str(sync.sync_id)},
            UpdateExpression="SET #state = :s, dataLocation = :dl, updatedAt = :ua",
            ExpressionAttributeValues={
                ":s": sync.state,
                ":dl": sync.data_location,
                ":ua": int(sync.updated_at.timestamp()),
            },
            ExpressionAttributeNames={
                "#state": "state"
            },  # state is a reserved keyword expression
        )

    def insert_or_update_items(self, sync: Sync) -> None:
        self._table.update_item(
            Key={"id": str(sync.sync_id)},
            UpdateExpression="SET #items = :its, #state = :s",
            ExpressionAttributeValues={
                ":its": [it.__dict__ for it in sync.items],
                ":s": sync.state,
            },
            ExpressionAttributeNames={
                "#items": "items",
                "#state": "state",
            },  # state & items are reserved keywords for expression
        )
