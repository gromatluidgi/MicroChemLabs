import pytest
from chembl_structure_pipeline import checker
from pandas import DataFrame
from rdkit import Chem
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


def test_mol2d_fixture(mol2d_frame: DataFrame):
    print(mol2d_frame.columns)
    mol: Mol = mol2d_frame["Molecule"][0]
    has_hygrogen: bool = MoleculeUtils.has_explicit_hygrogen(mol)

    assert has_hygrogen is False

    # Mol Block Extraction
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToMolBlock
    original_mol_block: str = Chem.MolToMolBlock(mol, True, -1, False)
    stero_mol_block: str = Chem.MolToMolBlock(mol, True, -1, True)
    kekulize_mol_block: str = Chem.MolToMolBlock(mol, True, -1, True)

    print("Has hygrogen: ", has_hygrogen)
    print("Original Mol Block: ", original_mol_block)
    print("Stereo Mol Block: ", stero_mol_block)
    print("Kekulize Mol Block: ", kekulize_mol_block)

    issues = checker.check_molblock(original_mol_block)
    print("Checking issues: ", issues)

    # Smiles
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToSmiles
    declared_smiles = mol2d_frame["PUBCHEM_OPENEYE_ISO_SMILES"][0]
    internal_can_smiles = Chem.MolToSmiles(mol)
    internal_can_no_stereo_smiles = Chem.MolToSmiles(mol, False)

    # Canonical smiles without stereo should not be affected for 2D molecule
    assert declared_smiles == internal_can_smiles == internal_can_no_stereo_smiles

    kekulized_smiles = Chem.MolToSmiles(mol, True, True, -1, False)
    kekulized_can_smiles = Chem.MolToSmiles(mol, True, True, -1, True)
    explicit_smiles = Chem.MolToSmiles(mol, True, True, -1, False, True, True)
    explicit_can_smiles = Chem.MolToSmiles(mol, True, True, -1, True, True, True)

    print("Declared SMILES: ", declared_smiles)
    print("Internal Canonical SMILES: ", internal_can_smiles)
    print("Internal Canonical SMILES without stereo: ", internal_can_no_stereo_smiles)
    print("Kekulized SMILES: ", kekulized_smiles)
    print("Kekulized Canonical SMILES: ", kekulized_can_smiles)
    print("Explicit SMILES: ", explicit_smiles)
    print("Explicit Canonical SMILES: ", explicit_can_smiles)


def test_mol3d_fixture(mol3d_frame: DataFrame):
    print(mol3d_frame.columns)
    mol: Mol = mol3d_frame["Molecule"][0]
    has_hygrogen: bool = MoleculeUtils.has_explicit_hygrogen(mol)

    # 3D mol should includes explicit hydrogen
    assert has_hygrogen is True

    # Mol Block Extraction
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToMolBlock
    original_mol_block: str = Chem.MolToMolBlock(mol, True, -1, False)
    stero_mol_block: str = Chem.MolToMolBlock(mol, True, -1, True)
    kekulize_mol_block: str = Chem.MolToMolBlock(mol, True, -1, True)

    print("Has hygrogen: ", has_hygrogen)
    print("Original Mol Block: ", original_mol_block)
    print("Stereo Mol Block: ", stero_mol_block)
    print("Kekulize Mol Block: ", kekulize_mol_block)

    issues = checker.check_molblock(original_mol_block)
    print("Checking issues: ", issues)

    # Smiles
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToSmiles
    declared_smiles = mol3d_frame["PUBCHEM_OPENEYE_ISO_SMILES"][0]
    internal_can_smiles = Chem.MolToSmiles(mol)
    internal_can_no_stereo_smiles = Chem.MolToSmiles(mol, False)

    # Canonical smiles without stereo should not be affected for 2D molecule
    # assert declared_smiles == internal_can_smiles == internal_can_no_stereo_smiles

    kekulized_smiles = Chem.MolToSmiles(mol, True, True, -1, False)
    kekulized_can_smiles = Chem.MolToSmiles(mol, True, True, -1, True)
    explicit_smiles = Chem.MolToSmiles(mol, True, True, -1, False, True, True)
    explicit_can_smiles = Chem.MolToSmiles(mol, True, True, -1, True, True, True)

    print("Declared SMILES: ", declared_smiles)
    print("Internal Canonical SMILES: ", internal_can_smiles)
    print("Internal Canonical SMILES without stereo: ", internal_can_no_stereo_smiles)
    print("Kekulized SMILES: ", kekulized_smiles)
    print("Kekulized Canonical SMILES: ", kekulized_can_smiles)
    print("Explicit SMILES: ", explicit_smiles)
    print("Explicit Canonical SMILES: ", explicit_can_smiles)


def test_molecules_from_sdf():
    sdf: str = PathUtils.join_path(
        PathUtils.get_parent_dir(__file__, 1), "_data", "metacetamol_3d.sdf"
    )

    result = []
    for mol in MoleculeUtils.molecules_from_sdf(sdf):
        result.append(mol)

    assert len(result) > 0
