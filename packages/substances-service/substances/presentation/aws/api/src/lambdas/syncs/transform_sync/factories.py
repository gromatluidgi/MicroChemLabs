import os

from shared.utils.path import PathUtils
from substances_core.domain.syncs.services.transformer import (
    SyncTransformer,
    SyncTransformerFactory,
)
from substances_etl.pubchem.transformer import (
    PubchemSubstanceTransformerOptions,
    PuchemSubstanceTransformer,
)


class SubstanceTransformerFactory(SyncTransformerFactory):
    """SubstanceTransformerFactory"""

    def from_provider(self, provider: str) -> SyncTransformer:
        if provider == "PUBCHEM":
            return PuchemSubstanceTransformer(
                PubchemSubstanceTransformerOptions(
                    data_dir=PathUtils.join_path(
                        PathUtils.root_path(), os.environ["EFS_MOUNT_DIR"], "syncs_data"
                    )
                )
            )
