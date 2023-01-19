from abc import ABC, abstractmethod
from typing import List, Optional

from shared.domain.repositories import Repository

from .objects import CanonicalSmiles, InchiKey
from .substance import Substance


class SubstanceReadRepository(Repository, ABC):
    @abstractmethod
    def find_by_inchi_key(self, inchi_key: InchiKey) -> Optional[Substance]:
        """
        Returns a substance with the given InchiKey.
        """
        raise NotImplementedError()

    @abstractmethod
    def find_by_smiles(self, smiles: CanonicalSmiles) -> Optional[Substance]:
        """
        Returns a substance with the given Canonical Smiles.
        """
        raise NotImplementedError()


class SubstanceWriteRepository(Repository, ABC):
    @abstractmethod
    def add(self, substance: Substance) -> None:
        """
        Add a new 'Substance' into the persistent store.
        """
        raise NotImplementedError()

    @abstractmethod
    def batch_insert_or_update(self, substances: List[Substance]) -> None:
        """
        Batch insert or update a batch of Subtance aggregates.
        """
        raise NotImplementedError()
