from _typeshed import Incomplete
from rdkit import Chem as Chem, Geometry as Geometry
from rdkit.Chem import rdDepictor as rdDepictor

def AlignDepict(mol, core, corePattern: Incomplete | None = ..., acceptFailure: bool = ...) -> None: ...
def initParser(): ...
def processArgs(args) -> None: ...
def main() -> None: ...