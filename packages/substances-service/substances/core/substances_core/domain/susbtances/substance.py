from dataclasses import dataclass, field
from datetime import datetime
from typing import List, Optional, Set

from shared.domain.aggregates import AggregateRoot

from .objects import (
    CanonicalSmiles,
    Formula,
    InchiKey,
    IupacName,
    PeptideSequence,
    SubstanceIdentifier,
    SubstanceName,
    SubstanceSynonym,
)


@dataclass
class Substance(AggregateRoot):
    """Substance aggregate root."""

    inchi_key: InchiKey

    formula: Formula
    """Standard chemical format (using Hill system) in protonated form to match reported charge."""

    mol_weight: float

    canonical_smiles: CanonicalSmiles

    fingerprints: List[str] = field(default_factory=list)

    identifiers: Set[SubstanceIdentifier] = field(default_factory=set)

    synonyms: Set[SubstanceSynonym] = field(default_factory=set)
    """List of alternative names of compound."""

    created_at: datetime = field(default_factory=datetime.utcnow)

    updated_at: Optional[datetime] = field(default=None)

    iupac_name: Optional[IupacName] = field(default=None)

    name: Optional[SubstanceName] = field(default=None)
    """Prefered name for the substance."""

    mass: Optional[float] = field(default=None)
    """Mass of compound or None when unknown."""

    avg_mol_weight: Optional[float] = field(default=None)
    """The average molecular weight of the compound."""

    mol2d: Optional[str] = field(default=None)
    """2D Molecule block"""

    mol3d: Optional[str] = field(default=None)
    """3D Molecule block."""

    sequence: Optional[PeptideSequence] = field(default=None)
    """Peptide sequence or None whne unknown."""

    standardized: bool = field(default=False)
    """
    Indicates wheiter this substance has been review by an expert or
    successfully pass an automated verification pipeline.
    """

    # def __post_init__(self):
    # self.add_domain_event(SubstanceCreatedEvent(self.inchi_key.key))

    def require_review(self) -> None:
        self.review_required = True
