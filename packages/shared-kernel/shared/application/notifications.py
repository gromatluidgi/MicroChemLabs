from abc import ABC, abstractmethod
from typing import Generic, TypeVar

from pydantic import BaseModel


class Notification(BaseModel):
    ...


TNotification = TypeVar("TNotification", None, Notification)


class NotificationHandler(ABC, Generic[TNotification]):
    """NotificationHandler"""

    @abstractmethod
    def handle(self, notification: TNotification) -> None:
        """Handle"""
        raise NotImplementedError


class NotificationPublisher(ABC, Generic[TNotification]):
    """NotificationPublisher"""

    @abstractmethod
    def notify(self, notification: TNotification) -> None:
        raise NotImplementedError
