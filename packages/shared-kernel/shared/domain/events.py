from abc import ABC
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import TypeVar


class EventScope(str, Enum):
    """EventScope"""

    INTERNAL = "internal"
    EXTERNAL = "external"
    BOTH = "both"


@dataclass
class DomainEvent:
    """DomainEvent"""

    aggregate_id: str = field(init=True)
    occured_on: datetime = field(default_factory=datetime.utcnow, init=False)
    scope: EventScope = field(default=EventScope.INTERNAL, init=False)


class EventStore(ABC):
    pass


# pylint: disable=invalid-name
TDomainEvent = TypeVar("TDomainEvent", None, DomainEvent)
"""TDomainEvent"""
