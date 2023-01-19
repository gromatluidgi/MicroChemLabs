import pytest
from pandas import DataFrame
from rdkit.Chem import Mol, PandasTools
from sciences.chemistry.molecules import MoleculeUtils
from shared.utils.path import PathUtils


@pytest.fixture(name="mol2d_frame")
def fixture_mol_2d():
    """Returns a pandas dataframe extracted from SDF file representing
    Metacetmol compound in 2D."""
    frame: DataFrame = PandasTools.LoadSDF(
        PathUtils.join_path(
            PathUtils.get_parent_dir(__file__, 1), "_data", "metacetamol_2d.sdf"
        ),
        idName="PUBCHEM_COMPOUND_CID",
        smilesName="PUBCHEM_OPENEYE_ISO_SMILES",
        molColName="Molecule",
    )
    return frame


@pytest.fixture(name="mol3d_frame")
def fixture_mol_3d():
    """Returns a pandas dataframe extracted from SDF file representing
    Metacetmol compound in 3D."""
    frame: DataFrame = PandasTools.LoadSDF(
        PathUtils.join_path(
            PathUtils.get_parent_dir(__file__, 1), "_data", "metacetamol_3d.sdf"
        ),
        idName="PUBCHEM_COMPOUND_CID",
        smilesName="PUBCHEM_OPENEYE_ISO_SMILES",
        molColName="Molecule",
    )
    return frame


@pytest.mark.parametrize(
    "frame, expected",
    [
        ("mol2d_frame", False),
        ("mol3d_frame", True),
    ],
)
def test_has_explicit_hydrogen(frame: str, expected: bool, request):
    mol: Mol = request.getfixturevalue(frame)["Molecule"][0]
    assert MoleculeUtils.has_explicit_hygrogen(mol) is expected


def test_get_canonical_smiles_from_mol():
    metamizole_3110 = MoleculeUtils.get_mol_from_smiles(
        "CC1=C(C(=O)N(N1C)C2=CC=CC=C2)N(C)CS(=O)(=O)[O-]"
    )
    metamizole_3111 = MoleculeUtils.get_mol_from_smiles(
        "CC1=C(C(=O)N(N1C)C2=CC=CC=C2)N(C)CS(=O)(=O)O"
    )

    stereo_can_smiles_3110 = MoleculeUtils.get_canonical_smiles_from_mol(
        metamizole_3110, True
    )
    can_smiles_3110 = MoleculeUtils.get_canonical_smiles_from_mol(
        metamizole_3110, False
    )

    stereo_can_smiles_3111 = MoleculeUtils.get_canonical_smiles_from_mol(
        metamizole_3111, True
    )
    can_smiles_3111 = MoleculeUtils.get_canonical_smiles_from_mol(
        metamizole_3111, False
    )

    assert stereo_can_smiles_3110 == can_smiles_3110
    assert stereo_can_smiles_3111 == can_smiles_3111
