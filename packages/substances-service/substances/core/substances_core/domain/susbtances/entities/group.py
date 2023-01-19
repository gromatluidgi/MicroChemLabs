"""_summary_
"""
from dataclasses import dataclass


@dataclass(frozen=True)
class FunctionalGroup:
    """SubstanceBehavior"""

    name: str
    """The name of the functional group."""

    _is_direct: bool = True

    def __post_init__(self) -> None:
        if self._is_direct:
            raise TypeError("Don't create instances of this class by hand!")

    @classmethod
    def amidine(cls):
        """
        Returns
        More info: https://en.wikipedia.org/wiki/Amidine
        """
        cls._is_direct = False
        cls.name = "Amidine"
        return cls

    @classmethod
    def carboxamide(cls):
        """
        Returns
        More info: https://en.wikipedia.org/wiki/Carboxamide
        """
        cls._is_direct = False
        cls.name = "Carboxamide"
        return cls
