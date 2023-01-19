import pytest
from sciences.chemistry.molecules import MoleculeUtils
from substances_importers.resolvers.cactus_resolver import CactusResolver


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(NC1=CC=CC(O)=C1)=O", "N-(3-Hydroxyphenyl)acetamide"),
        ("CC(=Nc1cccc(c1)O)O", None),
        ("CC(=O)NC1=CC(=CC=C1)O", "N-(3-Hydroxyphenyl)acetamide"),
        ("CC(=O)Nc1cccc(O)c1", "N-(3-Hydroxyphenyl)acetamide"),
    ],
)
def test_resolve_iupac_name(smiles, expected):
    iupac_name = CactusResolver.resolve_iupac_name(smiles)

    assert iupac_name == expected


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("CC(NC1=CC=CC(O)=C1)=O", "N-(3-Hydroxyphenyl)acetamide"),
        ("CC(=Nc1cccc(c1)O)O", "N-(3-Hydroxyphenyl)acetamide"),
        ("CC(=O)NC1=CC(=CC=C1)O", "N-(3-Hydroxyphenyl)acetamide"),
        ("CC(=O)Nc1cccc(O)c1", "N-(3-Hydroxyphenyl)acetamide"),
    ],
)
def test_resolve_iupac_name_can_smiles(smiles, expected):
    can_smiles = MoleculeUtils.get_canonical_smiles(smiles)
    iupac_name = CactusResolver.resolve_iupac_name(can_smiles)
    assert iupac_name == expected
