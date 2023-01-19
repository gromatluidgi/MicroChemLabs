from shared.application.events import DomainEventHandler
from shared.domain.aggregates import AggregateRoot
from shared.domain.events import DomainEvent
from shared.infrastructure.events import InternalDomainEventDispatcher


class FakeAggregate(AggregateRoot):
    ...


class FakeEvent(DomainEvent):
    name = "FakeEvent"


class FakeHandler(DomainEventHandler[FakeEvent]):
    def handle(self, event: FakeEvent) -> None:
        print("Fake Handle")


def test_register_handler():
    dispatcher = InternalDomainEventDispatcher()
    dispatcher.register_handler(FakeEvent, FakeHandler())


def test_dispatch():
    # Arrage
    dispatcher = InternalDomainEventDispatcher()
    dispatcher.register_handler(FakeEvent, FakeHandler())
    aggregate = FakeAggregate()
    aggregate.add_domain_event(FakeEvent("test"))

    # Act
    dispatcher.dispatch(aggregate)
