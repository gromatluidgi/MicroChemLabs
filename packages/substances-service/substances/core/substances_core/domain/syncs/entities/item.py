from dataclasses import dataclass, field
from typing import Optional, Union

from shared.domain.entities import Entity
from substances_core.domain.syncs.enums import SyncItemStatus


@dataclass
class SyncItem(Entity):
    """SyncItem"""

    token: str
    location: str
    name: Optional[str] = field(default=None)
    transformation: Optional[Union[str, bytes]] = field(default=None)
    status: SyncItemStatus = SyncItemStatus.PENDING

    def extracted_state(self) -> None:
        if self.status is not SyncItemStatus.PENDING:
            raise RuntimeError(
                "EXTRACTED status should be a successor of PENDING status."
            )
        self.status = SyncItemStatus.EXTRACTED

    def update_transformed_data(self, transformation: Union[str, bytes]) -> None:
        """
        Set the raw binary text represention of the transformed
        dataset.
        """
        if self.transformation is not None:
            raise RuntimeError("Cannot mutate transformation once realized.")

        self.transformation = transformation
        self.status = SyncItemStatus.TRANSFORMED
