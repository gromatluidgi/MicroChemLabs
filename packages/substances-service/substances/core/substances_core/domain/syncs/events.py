from dataclasses import dataclass

from shared.domain.events import DomainEvent


@dataclass
class SyncExecutedEvent(DomainEvent):
    """Event raised when a new sync is created."""

    provider: str


@dataclass
class SyncCheckedInEvent(DomainEvent):
    """SyncCheckedInEvent"""


@dataclass
class SyncCheckedOutEvent(DomainEvent):
    """SyncCheckedOutEvent"""
