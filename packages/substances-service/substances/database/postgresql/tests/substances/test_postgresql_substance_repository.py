import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from substances_core.domain.susbtances.objects import CanonicalSmiles, Formula, InchiKey
from substances_core.domain.susbtances.substance import Substance
from substances_postgresql.subtances.entities import SubstanceModel
from substances_postgresql.subtances.repository import (
    SubstancePostgresqlReadRepository,
    SubstancePostgresqlWriteRepository,
)


@pytest.fixture()
def database_engine():
    engine = create_engine(
        "postgresql://localhost/substances_test?user=admin&password=admin"
    )
    return engine


@pytest.fixture()
def database_session():
    engine = create_engine(
        "postgresql://localhost/substances_test?user=admin&password=admin"
    )
    session = sessionmaker(autocommit=False, autoflush=False, bind=engine)()
    try:
        yield session
    finally:
        session.close()


def test_count_row(database_session):
    count = database_session.query(SubstanceModel).count()
    print(count)


def test_add(database_session):
    substance = Substance(
        inchi_key=InchiKey("QLNWXBAGRTUKKI-UHFFFAOYSA-N"),
        formula=Formula("C8H9NO2"),
        mol_weight=0,
        canonical_smiles=CanonicalSmiles("COC(=O)C1=CC=CC=C1N"),
    )
    write_repo = SubstancePostgresqlWriteRepository(database_session)
    write_repo.add(substance)


def test_find_by_inchi_key(database_session):
    # Prepare
    inchikey = InchiKey("QLNWXBAGRTUKKI-UHFFFAOYSA-N")
    read_repo = SubstancePostgresqlReadRepository(database_session)

    # Act
    substance = read_repo.find_by_inchi_key(inchikey)

    assert substance is not None


def test_find_by_smiles(database_session):
    # Prepare
    smiles = CanonicalSmiles("COC(=O)C1=CC=CC=C1N")
    read_repo = SubstancePostgresqlReadRepository(database_session)

    # Act
    substance = read_repo.find_by_smiles(smiles)

    assert substance is not None
