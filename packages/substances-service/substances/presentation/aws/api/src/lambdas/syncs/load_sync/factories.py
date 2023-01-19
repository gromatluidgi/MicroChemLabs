import os

from shared.utils.path import PathUtils
from substances_core.domain.syncs.services.loader import SyncLoader, SyncLoaderFactory
from substances_etl.pubchem.loader import (
    PubchemSubstanceLoader,
    PubchemSubstanceLoaderOptions,
)


class SubstanceLoaderFactory(SyncLoaderFactory):
    """SubstanceLoaderFactory"""

    def from_provider(self, provider: str) -> SyncLoader:
        if provider == "PUBCHEM":
            return PubchemSubstanceLoader(
                PubchemSubstanceLoaderOptions(
                    data_dir=PathUtils.join_path(
                        PathUtils.root_path(), os.environ["EFS_MOUNT_DIR"], "syncs_data"
                    )
                )
            )
