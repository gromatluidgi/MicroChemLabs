from abc import ABC, abstractmethod
from typing import Callable, Generic, TypeVar

from pydantic import BaseModel
from shared.domain.events import DomainEvent, TDomainEvent


class Notification(BaseModel, Generic[TDomainEvent]):
    event: TDomainEvent


TNotification = TypeVar("TNotification", None, Notification[DomainEvent])


class NotificationHandler(ABC, Generic[TNotification]):
    @abstractmethod
    def handle(self, notification: TNotification) -> None:
        """Handle"""
        raise NotImplementedError


class NotificationPublisher(ABC, Generic[TNotification]):
    @abstractmethod
    def publish(self, notification: TNotification):
        raise NotImplementedError

    @abstractmethod
    def subscribe(self, event_name: str, func: Callable):
        raise NotImplementedError
