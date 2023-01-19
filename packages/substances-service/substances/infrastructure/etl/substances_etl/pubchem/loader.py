import gzip
import logging
import pickle  # nosec
from dataclasses import dataclass
from typing import List

import jsonpickle
from shared.utils.path import PathUtils  # type: ignore
from substances_core.domain.susbtances.substance import Substance
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.services.loader import SyncLoader

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PubchemSubstanceLoaderOptions:
    """PubchemSubstanceLoaderOptions"""

    data_dir: str


class PubchemSubstanceLoader(SyncLoader[Substance]):
    """PubchemSubstanceLoader"""

    def __init__(self, options: PubchemSubstanceLoaderOptions) -> None:
        self._options = options

    def load(self, item: SyncItem) -> List[Substance]:
        if item.transformation is None:
            raise ValueError

        if isinstance(item.transformation, str):
            pickle_path = PathUtils.join_path(
                self._options.data_dir, item.transformation
            )
            print(f"Load pickle: {pickle_path}")
            with open(pickle_path, "rb") as f_in:
                json = gzip.decompress(f_in.read())
                return jsonpickle.decode(pickle.loads(json))
        else:
            json = gzip.decompress(item.transformation)
            return jsonpickle.decode(pickle.loads(json))  # type: ignore #nosec
