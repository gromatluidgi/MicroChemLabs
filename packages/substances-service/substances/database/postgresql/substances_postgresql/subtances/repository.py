from datetime import datetime
from typing import List, Optional

from sqlalchemy.dialects.postgresql import insert
from sqlalchemy.orm import Session
from substances_core.domain.susbtances.objects import CanonicalSmiles, InchiKey
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.susbtances.substance import Substance
from substances_postgresql.subtances.entities import SubstanceModel


class SubstancePostgresqlReadRepository(SubstanceReadRepository):
    """
    SubstancePostgresqlReadRepository
    """

    def __init__(self, session: Session) -> None:
        self._session = session

    def find_by_inchi_key(self, inchi_key: InchiKey) -> Optional[Substance]:
        model: SubstanceModel = (
            self._session.query(SubstanceModel)
            .filter_by(_inchi_key=str(inchi_key))
            .one_or_none()
        )

        if model is None:
            return None

        return Substance(
            inchi_key=model.inchi_key,
            canonical_smiles=model.canonical_smiles,
            iupac_name=model.iupac_name,
            formula=model.formula,
            mol_weight=model.mol_weight,
            mass=model.mass,
            created_at=model.created_at,
            updated_at=model.updated_at,
        )

    def find_by_smiles(self, smiles: CanonicalSmiles) -> Optional[Substance]:
        model: SubstanceModel = (
            self._session.query(SubstanceModel)
            .filter_by(_canonical_smiles=str(smiles))
            .one_or_none()
        )

        if model is None:
            return None

        return Substance(
            inchi_key=model.inchi_key,
            canonical_smiles=model.canonical_smiles,
            iupac_name=model.iupac_name,
            formula=model.formula,
            mol_weight=model.mol_weight,
            mass=model.mass,
            created_at=model.created_at,
            updated_at=model.updated_at,
        )


class SubstancePostgresqlWriteRepository(SubstanceWriteRepository):
    """
    SubstancePostgresqlWriteRepository
    """

    __BATCH_MAX_CHUNK_SIZE: int = 1000

    def __init__(self, session: Session) -> None:
        self._session = session

    def add(self, substance: Substance) -> None:
        """
        Add a new 'Substance' into the persistent store.
        """
        substance_model = SubstanceModel.from_domain(substance=substance)
        self._session.add(substance_model)

    def batch_insert_or_update(self, substances: List[Substance]) -> None:
        for substance in substances:
            self._session.execute(
                insert(SubstanceModel)
                .values(SubstanceModel.dict_from_domain(substance))
                .on_conflict_do_update(
                    constraint=SubstanceModel.__table__.primary_key,
                    set_={"updated_at": datetime.utcnow()},
                )
            )
