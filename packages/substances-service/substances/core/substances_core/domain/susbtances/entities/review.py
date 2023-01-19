from dataclasses import field
from datetime import datetime
from typing import List


class Review:
    """Review"""

    id: str
    issues: List[str] = field(default_factory=list)
    created_at: datetime = field(default_factory=datetime.utcnow)
