from typing import List

from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.enums import SyncItemStatus


class SyncItemTranslator:
    @classmethod
    def from_record(cls, record) -> SyncItem:  # type: ignore
        return SyncItem(
            token=record["token"],
            location=record["location"],
            name=record["name"],
            transformation=record["transformation"]
            if record["transformation"] is not None
            else None,
            status=SyncItemStatus(record["status"]),
        )

    @classmethod
    def from_records(cls, records) -> List[SyncItem]:  # type: ignore
        syncs = []
        for record in records:
            syncs.append(SyncItemTranslator.from_record(record))
        return syncs
