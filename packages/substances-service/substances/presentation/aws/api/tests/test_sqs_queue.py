import boto3
from moto import mock_sqs
from shared.infrastructure.queue import Message, MessageType
from substances_aws_queue.queue import SQSMessageQueue


@mock_sqs
def test_queue_send_message():
    conn = boto3.resource("sqs", region_name="eu-west-3")
    conn.create_queue(QueueName="substances-syncs")

    queue = SQSMessageQueue("substances-syncs")
    queue.add_message(
        message=Message(id=1, type=MessageType.EVENT, body={"test": "test"})
    )
    queue.execute()
