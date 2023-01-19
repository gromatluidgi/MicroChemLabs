from typing import Any, Dict, Optional

from sqlalchemy import Column, Float, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from substances_core.domain.susbtances.objects import (
    CanonicalSmiles,
    Formula,
    InchiKey,
    IupacName,
)
from substances_core.domain.susbtances.substance import Substance
from substances_postgresql.mixins import TimestampMixin

Base = declarative_base()


class SubstanceModel(Base, TimestampMixin):
    """SubstanceModel"""

    __tablename__ = "substance"

    _inchi_key = Column("inchi_key", String, primary_key=True)
    _canonical_smiles = Column("canonical_smiles", String, index=True, nullable=False)
    _formula = Column("formula", String, index=True, nullable=False)
    _iupac_name = Column("iupac_name", String, index=True, nullable=True)
    mol_weight = Column(Float(asdecimal=True), nullable=False)
    mass = Column(String, nullable=True)
    mol2d = Column(String, nullable=True)
    mol3d = Column(String, nullable=True)

    @hybrid_property
    def inchi_key(self) -> InchiKey:
        return InchiKey(self._inchi_key)

    @inchi_key.setter
    def inchi_key(self, inchikey: InchiKey) -> None:
        if inchikey is not None:
            self._inchi_key = inchikey.key

    @hybrid_property
    def canonical_smiles(self) -> CanonicalSmiles:
        return CanonicalSmiles(self._canonical_smiles)

    @canonical_smiles.setter
    def canonical_smiles(self, smiles: CanonicalSmiles) -> None:
        if smiles is not None:
            self._canonical_smiles = smiles.value

    @hybrid_property
    def formula(self) -> Formula:
        return Formula(self._formula)

    @formula.setter
    def formula(self, formula: Formula) -> None:
        if formula is not None:
            self._formula = formula.molecular

    @hybrid_property
    def iupac_name(self) -> Optional[IupacName]:
        if self._iupac_name is not None:
            return IupacName(self._iupac_name)
        return None

    @iupac_name.setter
    def iupac_name(self, iupac_name: IupacName) -> None:
        if iupac_name is not None:
            self._iupac_name = iupac_name.value

    @classmethod
    def from_domain(cls, substance: Substance) -> "SubstanceModel":
        entity = cls()
        entity.inchi_key = substance.inchi_key
        entity.canonical_smiles = substance.canonical_smiles
        entity.formula = substance.formula
        entity.mol_weight = substance.mol_weight

        if substance.iupac_name is not None:
            entity.iupac_name = substance.iupac_name

        entity.mass = substance.mass
        entity.mol2d = substance.mol2d
        entity.mol3d = substance.mol3d
        return entity

    @classmethod
    def dict_from_domain(cls, substance: Substance) -> Dict[str, Any]:
        return {
            "inchi_key": substance.inchi_key.key,
            "canonical_smiles": substance.canonical_smiles.value,
            "mol_weight": substance.mol_weight,
            "iupac_name": substance.iupac_name.value
            if substance.iupac_name is not None
            else None,
            "formula": substance.formula.molecular,
            "mass": substance.mass,
            "mol2d": substance.mol2d,
            "mol3d": substance.mol3d,
        }
