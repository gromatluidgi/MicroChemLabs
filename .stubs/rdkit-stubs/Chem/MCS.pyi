from _typeshed import Incomplete

class MCSResult:
    numAtoms: Incomplete
    numBonds: Incomplete
    smarts: Incomplete
    completed: Incomplete
    def __init__(self, obj) -> None: ...
    def __nonzero__(self): ...

def FindMCS(mols, minNumAtoms: int = ..., maximize=..., atomCompare=..., bondCompare=..., matchValences=..., ringMatchesRingOnly: bool = ..., completeRingsOnly: bool = ..., timeout=..., threshold: Incomplete | None = ...): ...
