from abc import ABC
from dataclasses import dataclass
from typing import TypeVar


@dataclass
class Entity(ABC):
    """Entity"""


TEntity = TypeVar("TEntity", bound=Entity)  # pylint: disable=invalid-name
"""TEntity"""
