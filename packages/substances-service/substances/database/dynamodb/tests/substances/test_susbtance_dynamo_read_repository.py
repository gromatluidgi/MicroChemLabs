from moto import mock_dynamodb
from substances_core.domain.susbtances.objects import (
    CanonicalSmiles,
    Formula,
    InchiKey,
    IupacName,
)
from substances_core.domain.susbtances.substance import Substance
from substances_dynamodb.context import DynamoContext
from substances_dynamodb.substances.repository import (
    SubstanceDynamoReadRepository,
    SubstanceDynamoWriteRepository,
)


@mock_dynamodb
def test_batch_insert_or_update():
    context = DynamoContext()
    context._create_substances_table()
    read_repository = SubstanceDynamoReadRepository("substances")
    write_repository = SubstanceDynamoWriteRepository("substances")

    metacetamol = Substance(
        inchi_key=InchiKey("AOPRXJXHLWYPQR-UHFFFAOYSA-N"),
        formula=Formula("C8H9NO2"),
        canonical_smiles=CanonicalSmiles("NC(=O)COc1ccccc1"),
        mol_weight=151.165,
        iupac_name=IupacName("2-(phenoxy)acetamide"),
    )

    chloroquine = Substance(
        inchi_key=InchiKey("WHTVZRBIWZFKQO-UHFFFAOYSA-N"),
        formula=Formula("C18H26ClN3"),
        canonical_smiles=CanonicalSmiles("CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
        mol_weight=319.88,
        iupac_name=IupacName(
            "N'-(7-chloroquinolin-4-yl)-N,N-diethylpentane-1,4-diamine"
        ),
    )

    write_repository.batch_insert_or_update(substances=[metacetamol, chloroquine])

    assert read_repository.find_by_inchi_key(metacetamol.inchi_key) is not None
    assert read_repository.find_by_inchi_key(chloroquine.inchi_key) is not None
