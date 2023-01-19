from typing import Callable, Dict, Set, Type

from shared.application.events import DomainEventDispatcher, DomainEventHandler
from shared.domain.events import DomainEvent


class InternalDomainEventDispatcher(DomainEventDispatcher[DomainEvent]):
    """InternalDomainEventDispatcher"""

    def __init__(self) -> None:
        self._handlers: Dict[Type[DomainEvent], Set[Callable]] = {}

    def register_handler(
        self, event: Type[DomainEvent], handler: DomainEventHandler[DomainEvent]
    ) -> None:
        if event not in self._handlers:
            self._handlers[event] = set()

        self._handlers[event].add(handler)

    def dispatch(self, event: DomainEvent) -> None:
        self._do_dispatch(event)

    def _do_dispatch(self, event: DomainEvent) -> None:
        if type(event) in self._handlers:
            for handler in self._handlers[type(event)]:
                handler(event)
