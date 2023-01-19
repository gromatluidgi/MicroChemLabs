from abc import ABC, abstractmethod
from typing import List, Optional, Set

from substances_core.domain.syncs.entities.item import SyncItem


class SyncExtractor(ABC):
    """
    SyncExtractor
    """

    def get_data_location(self) -> str:
        return ""

    @abstractmethod
    def extract(self) -> List[SyncItem]:
        raise NotImplementedError


class SyncExtractorFactory(ABC):
    """SyncExtractingFactory"""

    # TODO: wrap params into config class
    @abstractmethod
    def from_provider(
        self,
        provider: str,
        data_location: Optional[str] = None,
        indexes: Set[str] = set(),
    ) -> SyncExtractor:
        raise NotImplementedError
