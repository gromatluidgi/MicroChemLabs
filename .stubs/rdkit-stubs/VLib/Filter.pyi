from _typeshed import Incomplete
from rdkit.VLib.Node import VLibNode as VLibNode

class FilterNode(VLibNode):
    def __init__(self, func: Incomplete | None = ..., negate: int = ..., **kwargs) -> None: ...
    def SetNegate(self, state) -> None: ...
    def Negate(self): ...
    def next(self): ...
