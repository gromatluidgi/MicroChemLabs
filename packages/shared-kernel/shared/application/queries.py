from abc import ABC, abstractmethod
from typing import Generic, List, TypeVar

from pydantic import BaseModel


class QueryBase(BaseModel):
    """QueryBase"""


class QueryResultBase(BaseModel):
    """QueryResultBase"""


T = TypeVar("T")


class QueryPaginationResult(QueryResultBase, Generic[T]):
    total: int
    items_per_page: int
    page: int
    items: List[T] = []


# pylint: disable=invalid-name
TQueryBase = TypeVar("TQueryBase", bound=QueryBase)
TQueryResultBase = TypeVar("TQueryResultBase", bound=QueryResultBase)


class QueryHandler(ABC, Generic[TQueryBase, TQueryResultBase]):
    """QueryHandler"""

    @abstractmethod
    def handle(self, query: TQueryBase) -> TQueryResultBase:
        """Handle"""
        raise NotImplementedError

    def __call__(self, query: TQueryBase) -> TQueryResultBase:
        """Make this class a callable."""
        return self.handle(query)
