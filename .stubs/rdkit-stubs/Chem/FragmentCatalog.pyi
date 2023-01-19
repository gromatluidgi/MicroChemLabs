from rdkit.Chem.rdfragcatalog import *
from _typeshed import Incomplete
from rdkit import Chem as Chem

def message(msg, dest=...) -> None: ...

class BitGainsInfo:
    id: int
    description: str
    gain: float
    nPerClass: Incomplete

def ProcessGainsFile(fileName, nToDo: int = ..., delim: str = ..., haveDescriptions: int = ...): ...
def BuildAdjacencyList(catalog, bits, limitInclusion: int = ..., orderLevels: int = ...): ...
def GetMolsMatchingBit(mols, bit, fps): ...
