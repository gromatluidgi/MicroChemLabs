import logging
import threading
from typing import List, Optional

import jsonpickle
import pika
from pydantic import BaseModel, Field
from shared.infrastructure.queue import Message, MessageQueue

LOGGER = logging.getLogger()


class RabbitMessageQueueOptions(BaseModel):
    """RabbitProducerOptions"""

    queue: str
    host: str
    port: int
    login: Optional[str] = Field(default=None)
    password: Optional[str] = Field(default=None)


class RabbitMessageQueue(MessageQueue):
    """RabbitMessageQueue"""

    def __init__(self, options: RabbitMessageQueueOptions) -> None:
        self._options = options
        self._queue: List[Message] = []
        self._connection = None
        self._channel = None

    @property
    def queue(self) -> List[Message]:
        return self._queue

    def add_message(self, message: Message):
        self._queue.append(message)

    def execute(self):
        connection_params = pika.ConnectionParameters(
            host=self._options.host,
            port=self._options.port,
        )
        self._connection = pika.SelectConnection(
            connection_params, on_open_callback=self._on_open_connection
        )
        self._open_connection()

    def _open_connection(self):
        """Perform async message publishing in a background thread (avoid main thread lock)"""
        try:
            iothread = threading.Thread(target=self._connection.ioloop.start, args=())
            iothread.start()
        except Exception as err:
            LOGGER.error(err)

    def _on_open_connection(self, _):
        self._channel = self._connection.channel(on_open_callback=self._on_open_channel)

    def _on_open_channel(self, _):
        self._channel.queue_declare(self._options.queue)
        for message in self._queue:
            self._channel.basic_publish(
                exchange="",
                routing_key=self._options.queue,
                body=jsonpickle.encode(message),
            )
            LOGGER.debug("Message sent: %r", message)
        self._close_connection()

    def _close_connection(self):
        self._connection.close()

    def _on_close_connection(self):
        pass
