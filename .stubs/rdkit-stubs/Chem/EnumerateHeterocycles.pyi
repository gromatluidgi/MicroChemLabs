from _typeshed import Incomplete
from collections.abc import Generator
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem

def GetHeterocycleReactionSmarts(): ...

REACTION_CACHE: Incomplete

def GetHeterocycleReactions(): ...
def EnumerateHeterocycles(inputmol, depth: Incomplete | None = ...) -> Generator[Incomplete, None, None]: ...
