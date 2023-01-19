from shared.infrastructure.queue import Message
from substances_core.domain.syncs.events import SyncExecutedEvent
from substances_rabbitmq.queue import RabbitMessageQueue, RabbitMessageQueueOptions


def test_send_empty_message():
    # Arrange
    message = Message(id=0, type="test")
    queue_options = RabbitMessageQueueOptions(
        queue="substances", host="localhost", port=5672
    )
    queue = RabbitMessageQueue(queue_options)
    queue.add_message(message)

    # Act
    queue.execute()


def test_send_message_with_body():
    # Arrange
    body = SyncExecutedEvent(aggregate_id="test", provider="test", state=1)
    message = Message(id=0, type="test", body=body)
    queue_options = RabbitMessageQueueOptions(
        queue="substances", host="localhost", port=5672
    )
    queue = RabbitMessageQueue(queue_options)
    queue.add_message(message)

    # Act
    queue.execute()
