from abc import ABC
from dataclasses import dataclass, field
from typing import List, TypeVar

from .entities import Entity
from .events import DomainEvent


@dataclass
class AggregateRoot(Entity, ABC):
    """AggregateRoot"""

    domain_events: List[DomainEvent] = field(default_factory=list, init=False)

    def add_domain_event(self, event: DomainEvent) -> None:
        self.domain_events.append(event)

    def clear_events(self) -> None:
        self.domain_events.clear()


# pylint: disable=invalid-name
TAggregateRoot = TypeVar("TAggregateRoot", bound=AggregateRoot)
"""TAggregateRoot"""
