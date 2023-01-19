from moto import mock_dynamodb
from substances_dynamodb.context import DynamoContext


@mock_dynamodb
def test_create_syncs_table():
    # Arrange
    context = DynamoContext()

    # Act
    table_desc = context._create_syncs_table()

    # Assert
    assert table_desc is not None


@mock_dynamodb
def test_create_substances_table():
    # Arrange
    context = DynamoContext()

    # Act
    table_desc = context._create_syncs_table()

    # Assert
    assert table_desc is not None
