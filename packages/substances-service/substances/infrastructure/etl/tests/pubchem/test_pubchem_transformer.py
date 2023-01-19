# pylint: disable=protected-access

import pytest
from pandas import DataFrame
from rdkit.Chem import PandasTools
from shared.utils.path import PathUtils
from substances_core.domain.syncs.entities.item import SyncItem
from substances_etl.pubchem.transformer import (
    PubchemSubstanceTransformerOptions,
    PuchemSubstanceTransformer,
)


@pytest.fixture(name="mol2d_frame")
def fixture_mol_2d():
    """Returns a pandas dataframe extracted from SDF file representing
    Metacetmol compound in 2 diemensions."""
    frame: DataFrame = PandasTools.LoadSDF(
        PathUtils.join_path(PathUtils.absolute_path(__file__), "metacetamol_2d.sdf"),
        idName="PUBCHEM_COMPOUND_CID",
        smilesName="PUBCHEM_OPENEYE_ISO_SMILES",
        molColName="Molecule",
    )
    return frame


@pytest.fixture(name="mol3d_frame")
def fixture_mol_3d():
    """Returns a pandas dataframe extracted from SDF file representing
    Metacetmol compound in 3 diemensions."""
    frame: DataFrame = PandasTools.LoadSDF(
        PathUtils.join_path(PathUtils.absolute_path(__file__), "metacetamol_3d.sdf"),
        idName="PUBCHEM_COMPOUND_CID",
        smilesName="PUBCHEM_OPENEYE_ISO_SMILES",
        molColName="Molecule",
    )
    return frame


@pytest.fixture
def transformer():
    data_dir = PathUtils.generate_temp_dir("MicroChemLabsTest")
    pubchem_transformer = PuchemSubstanceTransformer(
        PubchemSubstanceTransformerOptions(data_dir),
    )
    return pubchem_transformer


def test_transform(transformer: PuchemSubstanceTransformer):
    archive_path = PathUtils.join_path(
        PathUtils.generate_temp_dir("MicroChemLabsTest"),
        "Compound_127500001_128000000.sdf.gz",
    )
    item = SyncItem(
        token="363c446d401b0ce2fa0a0316faac39a5",
        location=archive_path,
        name="Compound_127500001_128000000.sdf.gz",
    )
    result = transformer.transform(item)
    assert result is not None
