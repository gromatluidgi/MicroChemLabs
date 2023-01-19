from abc import ABC, abstractmethod
from typing import Any, Generic

from shared.domain.events import TDomainEvent


class DomainEventHandler(Generic[TDomainEvent], ABC):
    """DomainEventHandler"""

    @abstractmethod
    def handle(self, event: TDomainEvent) -> None:
        """Handle"""
        raise NotImplementedError

    def __call__(self, event: TDomainEvent) -> Any:
        """Make this class a callable."""
        return self.handle(event)


class DomainEventDispatcher(Generic[TDomainEvent], ABC):
    """DomainEventDispatcher"""

    @abstractmethod
    def dispatch(self, event: TDomainEvent) -> None:
        """Dispatch the events occured for a specific aggregate."""
        raise NotImplementedError
