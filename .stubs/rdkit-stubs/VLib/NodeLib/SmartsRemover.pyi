from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.VLib.Transform import TransformNode as TransformNode

class SmartsRemover(TransformNode):
    def __init__(self, patterns=..., wholeFragments: int = ..., **kwargs) -> None: ...
    def transform(self, cmpd): ...

biggerTest: str
__test__: Incomplete
