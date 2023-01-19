import logging
from typing import Optional

import requests
from substances_core.domain.susbtances.objects import IupacName

LOGGER = logging.getLogger(__name__)


class CactusResolver:

    __CACTUS_IUPAC_NAME_URI: str = (
        "https://cactus.nci.nih.gov/chemical/structure/{0}/iupac_name"
    )

    @staticmethod
    def resolve_iupac_name(smiles: str) -> Optional[IupacName]:
        url = CactusResolver.__CACTUS_IUPAC_NAME_URI.format(smiles)
        response = requests.get(url)
        try:
            response.raise_for_status()
        except Exception as err:
            LOGGER.error(err)
            return None
        return IupacName(response.text)
