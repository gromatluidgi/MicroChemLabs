import pytest
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import Mol, PandasTools
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


def test_mol2d_to_inchi(mol2d_frame: str):
    mol: Mol = mol2d_frame["Molecule"][0]

    extracted_inchi = mol2d_frame["PUBCHEM_IUPAC_INCHI"][0]
    extracted_inchi_key = mol2d_frame["PUBCHEM_IUPAC_INCHIKEY"][
        0
    ]  # 'QLNWXBAGRTUKKI-UHFFFAOYSA-N'

    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.MolToInchiKey(mol)
    print("Inchi: ", inchi)
    print("InchiKey: ", inchikey)

    assert extracted_inchi == inchi
    assert extracted_inchi_key == inchikey


def test_mol3d_to_inchi(mol3d_frame: str):
    # Arrange
    mol: Mol = mol3d_frame["Molecule"][0]

    # Act
    inchi = Chem.MolToInchi(mol)
    inchikey = Chem.MolToInchiKey(mol)
    print("Inchi: ", inchi)
    print("InchiKey: ", inchikey)

    # Assert
    assert inchi == "InChI=1S/C8H9NO2/c1-6(10)9-7-3-2-4-8(11)5-7/h2-5,11H,1H3,(H,9,10)"
    assert inchikey == "QLNWXBAGRTUKKI-UHFFFAOYSA-N"
