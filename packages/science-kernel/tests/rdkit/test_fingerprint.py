import pytest
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import Mol, PandasTools
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.DataStructs.cDataStructs import (
    BitVectToBinaryText,
    BitVectToFPSText,
    BitVectToText,
    ConvertToNumpyArray,
)
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


def test_path_fingerprint(mol2d_frame):
    mol: Mol = mol2d_frame["Molecule"][0]
    fingerprint = Chem.RDKFingerprint(mol)
    print("Binary Text Fingerprint: ", BitVectToBinaryText(fingerprint))
    print("FPS Text Fingerprint: ", BitVectToFPSText(fingerprint))
    print("Text Fingerprint: ", BitVectToText(fingerprint))


def test_morgan_fingerprint():
    ...
