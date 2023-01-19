from abc import ABC, abstractmethod
from enum import Enum
from typing import Any, List

from substances_core.domain.susbtances.substance import Substance


class ImporterDataType(str, Enum):
    FILE = "file"
    STRING = "string"
    BINARY = "binary"


class SubstanceImporterReader(ABC):
    @abstractmethod
    def read(self, data: Any) -> List[Substance]:
        raise NotImplementedError


class SubstanceImporter(ABC):
    @abstractmethod
    def load(self, data: Any) -> List[Substance]:
        raise NotImplementedError


class SubstanceImporterConfiguration:
    def __init__(self) -> None:
        pass


class SubstanceImporterFactory(ABC):
    @abstractmethod
    def importer_by_format(
        self, source_type: ImporterDataType, format: str
    ) -> SubstanceImporter:
        raise NotImplementedError
