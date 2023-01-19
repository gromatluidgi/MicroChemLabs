import logging
from typing import List

import boto3
import jsonpickle
from shared.infrastructure.queue import Message, MessageQueue

LOGGER = logging.getLogger()


class SQSMessageQueue(MessageQueue):
    """RabbitMessageQueue"""

    def __init__(self, queue_name: str) -> None:
        self._messages: List[Message] = []
        self._queue_service = boto3.resource("sqs")
        self._queue = self._queue_service.Queue(queue_name)

    @property
    def queue(self) -> List[Message]:
        return self._messages

    def add_message(self, message: Message):
        self._messages.append(message)

    def execute(self):
        for message in self._messages:
            response = self._queue.send_message(MessageBody=jsonpickle.encode(message))
            LOGGER.debug("Message sent: %r", message)
            LOGGER.debug("Message response: %r", response)
