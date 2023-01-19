from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Graphs as Graphs, rdMolDescriptors as rdMolDescriptors, rdchem as rdchem
from rdkit.ML.InfoTheory import entropy as entropy

ptable: Incomplete
hallKierAlphas: Incomplete

def Ipc(mol, avg: int = ..., dMat: Incomplete | None = ..., forceDMat: int = ...): ...

HallKierAlpha: Incomplete
Kappa1: Incomplete
Kappa2: Incomplete
Kappa3: Incomplete

def Chi0(mol): ...
def Chi1(mol): ...

Chi0v: Incomplete
Chi1v: Incomplete
Chi2v: Incomplete
Chi3v: Incomplete
Chi4v: Incomplete
ChiNv_: Incomplete
Chi0n: Incomplete
Chi1n: Incomplete
Chi2n: Incomplete
Chi3n: Incomplete
Chi4n: Incomplete
ChiNn_: Incomplete

def BalabanJ(mol, dMat: Incomplete | None = ..., forceDMat: int = ...): ...
def BertzCT(mol, cutoff: int = ..., dMat: Incomplete | None = ..., forceDMat: int = ...): ...
