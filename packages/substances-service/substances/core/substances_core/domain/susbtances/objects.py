import re
from dataclasses import dataclass

_INCHIKEY_REGEX = r"^([0-9A-Z\-]+)$"
_INCHIKEY_PATTERN = re.compile(_INCHIKEY_REGEX)
_CASNUMBER_REGEX = r""
_CASNUMBER_PATTERN = re.compile(_CASNUMBER_REGEX)


@dataclass(frozen=True, eq=True)
class InchiKey:
    """
    InChIKey is a fixed-length format directly derived from InChI. It is based on a strong hash (SHA-256 algorithm)
    of an InChI string.

    The 27 characters long InChIKey is made of three parts connected by hyphens.
    The first part is 14 characters long and is based on the connectivity and proton layers of an InChI string.
    The second part, contains 9 characters that are related to all other InChI layers (isotopes, stereochemistry, etc.)
    and also contains the version of InChI and its standard/non-standard property in the last two characters.
    The third part is one letter, describing the (de)protonation layer of the original InChI.

    http://inchi.info/inchikey_overview_en.html
    """

    key: str

    def __post_init__(self) -> None:
        # Basic check before trying regex match
        if (len(self.key) != 27 or self.key[14] != "-") or _INCHIKEY_PATTERN.match(
            self.key
        ) is None:
            raise ValueError("InchiKey must respect a valid standard format.")

    def __str__(self) -> str:
        return self.key


class Inchi:
    value: str


@dataclass(frozen=True, eq=True)
class CasNumber:
    value: str


@dataclass(frozen=True, eq=True)
class CanonicalSmiles:
    value: str

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, eq=True)
class SubstanceIdentifier:
    scheme: str
    value: str


@dataclass(frozen=True, eq=True)
class SubstanceName:
    value: str

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, eq=True)
class IupacName:
    value: str

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, eq=True)
class SubstanceSynonym:
    value: str


@dataclass(frozen=True)
class PeptideSequence:
    value: str

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True)
class Formula:
    """Formula represent the chemical composition of a substance."""

    molecular: str

    def __str__(self) -> str:
        return self.molecular


@dataclass(frozen=True, eq=True)
class Issue:
    """Issue"""

    code: str
    description: str


@dataclass(frozen=True, eq=True)
class Fingerprint:
    pass


@dataclass(frozen=True, eq=True)
class Vendor:
    name: str
    product_id: str
