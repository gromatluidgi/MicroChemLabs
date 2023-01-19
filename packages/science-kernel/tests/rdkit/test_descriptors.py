import pytest
from pandas import DataFrame
from rdkit.Chem import Descriptors, Mol, PandasTools
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
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


# https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html#rdkit.Chem.rdMolDescriptors.CalcMolFormula
@pytest.mark.parametrize(
    "frame_fixture",
    [
        ("mol2d_frame"),
        ("mol3d_frame"),
    ],
)
def test_calc_mol_formula(frame_fixture: str, request):
    frame = request.getfixturevalue(frame_fixture)
    mol: Mol = frame["Molecule"][0]

    formula = CalcMolFormula(mol)

    assert formula == "C8H9NO2"


@pytest.mark.parametrize(
    "frame_fixture",
    [
        ("mol2d_frame"),
        ("mol3d_frame"),
    ],
)
def test_get_mol_weigth(frame_fixture: str, request):
    frame = request.getfixturevalue(frame_fixture)
    mol: Mol = frame["Molecule"][0]

    weight = Descriptors.MolWt(mol)

    assert round(weight, 2) == 151.16
