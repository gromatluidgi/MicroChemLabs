# See:https://github.com/pika/pika/blob/main/examples/asynchronous_consumer_example.py
import logging
import threading
from typing import Any

import jsonpickle
import pika
from shared.infrastructure.queue import Message, MessageType

from ..config import settings
from .dependencies import internal_command_dispatcher


class RabbitConsumer:
    """RabbitConsumer"""

    log = logging.getLogger(__name__)
    connection = None
    channel = None

    @classmethod
    def run(cls):
        """Start rabbitmq consuming"""
        cls._connect()

    @classmethod
    def _connect(cls):
        cls.connection = pika.SelectConnection(
            pika.ConnectionParameters(host=settings.RABBITMQ_HOST),
            on_open_callback=cls._on_connection_open,
        )
        # Execute the ioloop into a background process
        iothread = threading.Thread(target=cls.connection.ioloop.start, args=())
        iothread.start()

    @classmethod
    def _on_connection_open(cls, _unused_connection):
        """This method is called by pika once the connection to RabbitMQ has
        been established. It passes the handle to the connection object in
        case we need it, but in this case, we'll just mark it unused.
        :param pika.SelectConnection _unused_connection: The connection
        """
        cls.log.info("Connection opened")
        cls._open_channel()

    @classmethod
    def _open_channel(cls):
        cls.connection.channel(on_open_callback=cls._on_channel_open)

    @classmethod
    def _on_channel_open(cls, channel):
        """This method is invoked by pika when the channel has been opened.
        The channel object is passed in so we can make use of it.
        Since the channel is now open, we'll declare the exchange to use.
        :param pika.channel.Channel channel: The channel object
        """
        cls.log.info("Channel opened")
        cls.channel = channel
        cls._setup_queue()

    @classmethod
    def _setup_queue(cls):
        cls.channel.queue_declare(
            queue=settings.RABBITMQ_NAME, callback=cls._on_queue_declareok
        )

    @classmethod
    def _on_queue_declareok(cls, _):
        """Method invoked by pika when the Queue.Declare RPC call made in
        setup_queue has completed. In this method we will bind the queue
        and exchange together with the routing key by issuing the Queue.Bind
        RPC command. When this command is complete, the on_bindok method will
        be invoked by pika.
        :param pika.frame.Method _unused_frame: The Queue.DeclareOk frame
        :param str|unicode userdata: Extra user data (queue name)
        """
        cls._start_consuming()

    @classmethod
    def _start_consuming(cls):
        cls.channel.basic_consume(
            queue=settings.RABBITMQ_NAME,
            on_message_callback=RabbitConsumer._message_callback,
            auto_ack=True,
        )

    @classmethod
    def _message_callback(cls, _: Any, __: Any, ___: Any, body: Any):
        cls.log.info("Received message raw body: %r", body)
        try:
            message: Message = jsonpickle.decode(body)
            if message.type == MessageType.COMMAND:
                iothread = threading.Thread(
                    target=internal_command_dispatcher().dispatch, args=(message.body,)
                )
                iothread.start()
            cls.log.info("Deserialized Message: %r", message)
        except Exception as err:
            cls.log.error(err)
