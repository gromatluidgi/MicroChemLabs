# pylint: disable=protected-access
import pytest
from shared.utils.path import PathUtils
from substances_etl.pubchem.extractor import (
    PubchemExtractorOptions,
    PuchemSubstanceExtractor,
)


@pytest.fixture
def extractor():
    host = "ftp.ncbi.nlm.nih.gov"
    data_location = "pubchem/Compound/CURRENT-Full/SDF"
    data_dir = PathUtils.generate_temp_dir("MicroChemLabsTest")
    pubchem_extractor = PuchemSubstanceExtractor(
        PubchemExtractorOptions(
            host,
            data_location,
            data_dir,
            indexes=["Compound_127500001_128000000.sdf.gz"],
        ),
    )
    return pubchem_extractor


def test_extract(extractor: PuchemSubstanceExtractor):
    result = extractor.extract()
    assert result is not None
