from datetime import datetime, timedelta
from typing import Optional, Set

from shared.utils.path import PathUtils
from substances_core.domain.syncs.services.extractor import (
    SyncExtractor,
    SyncExtractorFactory,
)
from substances_core.domain.syncs.services.loader import SyncLoader, SyncLoaderFactory
from substances_core.domain.syncs.services.transformer import (
    SyncTransformer,
    SyncTransformerFactory,
)
from substances_etl.pubchem.extractor import (
    PubchemExtractorOptions,
    PuchemSubstanceExtractor,
)
from substances_etl.pubchem.loader import (
    PubchemSubstanceLoader,
    PubchemSubstanceLoaderOptions,
)
from substances_etl.pubchem.transformer import (
    PubchemSubstanceTransformerOptions,
    PuchemSubstanceTransformer,
)

from .config import settings


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

            # TODO: use settings || env
            return PuchemSubstanceExtractor(
                PubchemExtractorOptions(
                    host=settings.PUBCHEM_FTP_HOST,
                    data_location=data_location,
                    data_dir=PathUtils.join_path(
                        PathUtils.root_path(), "_data", "MicroChemLabs"
                    ),
                    items_limit=1000,
                ),
            )


class SubstanceTransformerFactory(SyncTransformerFactory):
    """SubstanceTransformerFactory"""

    # TODO: use settings || env
    def from_provider(self, provider: str) -> SyncTransformer:
        if provider == "PUBCHEM":
            return PuchemSubstanceTransformer(
                PubchemSubstanceTransformerOptions(
                    data_dir=PathUtils.join_path(
                        PathUtils.root_path(), "_data", "MicroChemLabs"
                    )
                )
            )


class SubstanceLoaderFactory(SyncLoaderFactory):
    """SubstanceLoaderFactory"""

    # TODO: use settings || env
    def from_provider(self, provider: str) -> SyncLoader:
        if provider == "PUBCHEM":
            return PubchemSubstanceLoader(
                PubchemSubstanceLoaderOptions(
                    data_dir=PathUtils.join_path(
                        PathUtils.root_path(), "_data", "MicroChemLabs"
                    )
                )
            )
