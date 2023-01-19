from enum import Enum


class SyncItemStatus(str, Enum):
    """SyncItemStatus"""

    PENDING = "pending"
    """Default status for an item awaiting extraction."""

    EXTRACTED = "extracted"
    """Indicates that item is fetched and extracted into data storage."""

    TRANSFORMED = "transformed"
    """Indicates that extracted item has been transformed and ready for load."""

    LOADED = "loaded"

    ERROR = "error"


class SyncState(str, Enum):
    """SyncState"""

    PENDING = "pending"
    CHECKIN = "checkin"
    CHECKOUT = "checkout"
    COMPLETED = "completed"
    ERROR = "error"


class SyncScope(str, Enum):
    """SyncScope"""

    SUBSTANCE = "substance"
    VENDOR = "vendor"
