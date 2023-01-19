from _typeshed import Incomplete
from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Chem.rdmolfiles import SDMolSupplier as SDMolSupplier, SmilesMolSupplier as SmilesMolSupplier

class InputFormat:
    SMARTS: str
    MOL: str
    SMILES: str

class SaltRemover:
    defnFilename: Incomplete
    defnData: Incomplete
    salts: Incomplete
    defnFormat: Incomplete
    def __init__(self, defnFilename: Incomplete | None = ..., defnData: Incomplete | None = ..., defnFormat=...) -> None: ...
    def StripMol(self, mol, dontRemoveEverything: bool = ..., sanitize: bool = ...): ...
    def StripMolWithDeleted(self, mol, dontRemoveEverything: bool = ...): ...
    def __call__(self, mol, dontRemoveEverything: bool = ...): ...
