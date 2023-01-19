from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Crippen as Crippen, rdMolDescriptors as rdMolDescriptors, rdPartialCharges as rdPartialCharges

ptable: Incomplete
bondScaleFacts: Incomplete
mrBins: Incomplete

def pySMR_VSA_(mol, bins: Incomplete | None = ..., force: int = ...): ...

SMR_VSA_: Incomplete
logpBins: Incomplete

def pySlogP_VSA_(mol, bins: Incomplete | None = ..., force: int = ...): ...

SlogP_VSA_: Incomplete
chgBins: Incomplete

def pyPEOE_VSA_(mol, bins: Incomplete | None = ..., force: int = ...): ...

PEOE_VSA_: Incomplete

def pyLabuteASA(mol, includeHs: int = ...): ...

LabuteASA: Incomplete
TPSA: Incomplete
