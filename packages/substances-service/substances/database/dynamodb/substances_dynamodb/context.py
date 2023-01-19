import boto3


class DynamoContext:
    def __init__(self) -> None:
        pass

    def _create_syncs_table(self) -> None:
        client = boto3.client("dynamodb", region_name="eu-west-3")
        client.create_table(
            TableName="syncs",
            KeySchema=[
                {"AttributeName": "id", "KeyType": "HASH"},  # Partition Key
            ],
            AttributeDefinitions=[
                {"AttributeName": "id", "AttributeType": "S"},
                {"AttributeName": "provider", "AttributeType": "S"},
                {"AttributeName": "createdAt", "AttributeType": "N"},
            ],
            GlobalSecondaryIndexes=[
                {
                    "IndexName": "SyncProviderIndex",
                    "KeySchema": [
                        {"AttributeName": "provider", "KeyType": "HASH"},
                        {"AttributeName": "createdAt", "KeyType": "RANGE"},  # Sort key
                    ],
                    "Projection": {
                        "ProjectionType": "INCLUDE",
                        "NonKeyAttributes": ["state", "dataLocation"],
                    },
                    "ProvisionedThroughput": {
                        "ReadCapacityUnits": 1,
                        "WriteCapacityUnits": 1,
                    },
                },
            ],
            ProvisionedThroughput={"ReadCapacityUnits": 1, "WriteCapacityUnits": 1},
        )
        return client.describe_table(TableName="syncs")["Table"]

    def _create_substances_table(self) -> None:
        client = boto3.client("dynamodb", region_name="eu-west-3")
        client.create_table(
            TableName="substances",
            KeySchema=[
                {"AttributeName": "id", "KeyType": "HASH"},
            ],
            AttributeDefinitions=[
                {"AttributeName": "id", "AttributeType": "S"},
            ],
            ProvisionedThroughput={"ReadCapacityUnits": 1, "WriteCapacityUnits": 1},
        )
        return client.describe_table(TableName="substances")["Table"]
