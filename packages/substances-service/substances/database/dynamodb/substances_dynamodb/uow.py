from typing import Callable, Dict, List

from shared.domain.events import DomainEvent
from shared.domain.repositories import Repository
from shared.domain.uow import UnitOfWork
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.syncs.repository import (
    SyncReadRepository,
    SyncWriteRepository,
)
from substances_dynamodb.substances.repository import (
    SubstanceDynamoReadRepository,
    SubstanceDynamoWriteRepository,
)
from substances_dynamodb.syncs.repository import (
    SyncDynamoReadRepository,
    SyncDynamoWriteRepository,
)


class DynamoUnitOfWork(UnitOfWork):
    """Unit of work implementation for DynamoDB."""

    def __init__(self, table_names: Dict[type, str]) -> None:
        self._table_names = table_names
        self._repositories: Dict[type, Callable[[str], Repository]] = {
            SyncReadRepository: SyncDynamoReadRepository,
            SyncWriteRepository: SyncDynamoWriteRepository,
            SubstanceReadRepository: SubstanceDynamoReadRepository,
            SubstanceWriteRepository: SubstanceDynamoWriteRepository,
        }

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    @property
    def repositories(self) -> Dict[str, Callable[[], Repository]]:
        return self._repositories

    def get_repository(self, repository: type[Repository]) -> Repository:
        if repository not in self.repositories:
            raise Exception("No repository found!")

        if repository not in self._table_names:
            raise Exception("No name found for table.")

        table_name = self._table_names[repository]

        return self._repositories[repository](table_name)

    def commit(self, events: List[DomainEvent] = []):
        """Transaction are handled directly within repository"""
        pass
