import uuid
from abc import ABC, abstractmethod
from datetime import datetime
from enum import Enum
from typing import Any, List, Optional, Union
from uuid import UUID

from pydantic import BaseModel, Json
from shared.application.commands import CommandBase


class MessageType(str, Enum):
    COMMAND = "command"
    EVENT = "event"
    RAW = "raw"


class Message(BaseModel):
    """Message"""

    id: Union[int, str, UUID] = uuid.uuid4()
    type: MessageType
    created_at: datetime = datetime.utcnow()
    expires_at: Optional[datetime] = None
    body: Union[CommandBase, Json[Any], None] = None

    class Config:
        allow_mutation = False


class MessageQueue(ABC):
    @property
    @abstractmethod
    def queue(self) -> List[Message]:
        raise NotImplementedError

    @abstractmethod
    def add_message(self, message: Message) -> None:
        """Add a message to the queue."""
        raise NotImplementedError

    @abstractmethod
    def execute(self) -> None:
        """Process all messages in the queue and return count."""
        raise NotImplementedError
