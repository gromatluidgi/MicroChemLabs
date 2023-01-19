import calendar
import logging
from datetime import datetime
from typing import List, Optional

import boto3
from botocore.exceptions import ClientError
from substances_core.domain.susbtances.objects import (
    CanonicalSmiles,
    Formula,
    InchiKey,
    IupacName,
)
from substances_core.domain.susbtances.repository import (
    SubstanceReadRepository,
    SubstanceWriteRepository,
)
from substances_core.domain.susbtances.substance import Substance

LOGGER = logging.getLogger(__name__)


class SubstanceDynamoReadRepository(SubstanceReadRepository):
    """
    SubstanceDynamoReadRepository
    """

    def __init__(self, table_name: str) -> None:
        self._dyn_resource = boto3.resource("dynamodb")
        self._table = self._dyn_resource.Table(table_name)

    def find_by_inchi_key(self, inchi_key: InchiKey) -> Optional[Substance]:
        """
        Gets sync data from the data store for a specific sync.
        """
        try:
            response = self._table.get_item(Key={"id": inchi_key.key})
        except ClientError as err:
            LOGGER.error(
                "Couldn't get sync %s from table %s. Here's why: %s: %s",
                str(inchi_key),
                self._table.name,
                err.response["Error"]["Code"],
                err.response["Error"]["Message"],
            )
            return None
        else:
            item = response["Item"]
            return Substance(
                inchi_key=InchiKey(item["id"]),
                formula=Formula(item["formula"]),
                mol_weight=float(item["molWeight"]),
                canonical_smiles=CanonicalSmiles(item["canonicalSmiles"]),
                iupac_name=IupacName(item["iupacName"]) if item["iupacName"] else None,
                created_at=datetime.fromtimestamp(item["createdAt"]),
            )

    def find_by_smiles(self, smiles: CanonicalSmiles) -> Optional[Substance]:
        return super().find_by_smiles(smiles)


class SubstanceDynamoWriteRepository(SubstanceWriteRepository):
    """
    SyncDynamoWriteRepository
    """

    def __init__(self, table_name: str) -> None:
        self._table_name = table_name
        self._dyn_resource = boto3.resource("dynamodb")
        self._table = self._dyn_resource.Table(table_name)

    def add(self, substance: Substance) -> None:
        """
        Adds a substance to the DynamoDB table.
        """
        self._table.put_item(
            Item={
                "id": str(substance.inchi_key),
                "iupacName": str(substance.iupac_name)
                if substance.iupac_name
                else None,
                "formula": substance.formula.molecular,
                "molWeight": substance.mol_weight,
                "canonicalSmiles": substance.canonical_smiles.value,
                "createdAt": calendar.timegm(substance.created_at.utctimetuple()),
            }
        )

    def __chunk(self, substances: List[Substance], chunk_size):
        for i in range(0, len(substances), chunk_size):
            yield substances[i : (i + chunk_size)]

    def __batch_item_structure(self, substance: Substance):
        return {
            "id": str(substance.inchi_key),
            "iupacName": str(substance.iupac_name) if substance.iupac_name else None,
            "formula": substance.formula.molecular,
            "molWeight": str(substance.mol_weight),
            "canonicalSmiles": substance.canonical_smiles.value,
            "createdAt": calendar.timegm(substance.created_at.utctimetuple()),
        }

    def batch_insert_or_update(self, substances: List[Substance]) -> None:
        """
        Batch insert or update Substance aggregates.
        """
        # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/dynamodb.html#DynamoDB.Client.batch_write_item

        count = 0
        for chunk in self.__chunk(substances, 25):
            items = {
                self._table_name: [
                    {"PutRequest": {"Item": self.__batch_item_structure(item)}}
                    for item in chunk
                ]
            }
            response = self._dyn_resource.batch_write_item(RequestItems=items)
            print(f"Unprocessed items: {response['UnprocessedItems']}")
            print(f"Chunk {count}")
            count += 1
