import os
from datetime import datetime, timedelta
from typing import Optional, Set

from shared.utils.path import PathUtils
from substances_core.domain.syncs.services.extractor import (
    SyncExtractor,
    SyncExtractorFactory,
)
from substances_etl.pubchem.extractor import (
    PubchemExtractorOptions,
    PuchemSubstanceExtractor,
)

PUBCHEM_FTP_HOST: str = "ftp.ncbi.nlm.nih.gov"  # TODO: get from env


class SubstanceExtractorFactory(SyncExtractorFactory):
    """SubstanceExtractorFactory"""

    def from_provider(
        self,
        provider: str,
        data_location: Optional[str] = None,
        indexes: Set[str] = set(),
    ) -> SyncExtractor:
        if provider == "PUBCHEM":
            if data_location is None:
                yesterday = datetime.now() - timedelta(days=1)
                data_location = (
                    f'pubchem/Compound/Daily/{yesterday.strftime("%Y-%m-%d")}/SDF'
                )

            return PuchemSubstanceExtractor(
                PubchemExtractorOptions(
                    PUBCHEM_FTP_HOST,
                    data_location,
                    PathUtils.join_path(
                        PathUtils.root_path(), os.environ["EFS_MOUNT_DIR"], "syncs_data"
                    ),
                    indexes=indexes,
                ),
            )
        raise ValueError("Unknown provider")
