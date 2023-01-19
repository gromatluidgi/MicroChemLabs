from _typeshed import Incomplete
from rdkit import Chem as Chem, DataStructs as DataStructs
from rdkit.Chem import rdqueries as rdqueries

rdkitFpParams: Incomplete
FTYPE_ACYCLIC: str
FTYPE_CYCLIC: str
FTYPE_CYCLIC_ACYCLIC: str
ACYC_SMARTS: Incomplete
CYC_SMARTS: Incomplete
cSma1: Incomplete
cSma2: Incomplete
dummyAtomQuery: Incomplete

def delete_bonds(mol, bonds, ftype, hac): ...
def select_fragments(fragments, ftype, hac): ...
def isValidRingCut(mol): ...
def generate_fraggle_fragmentation(mol, verbose: bool = ...): ...
def atomContrib(subs, mol, tverskyThresh: float = ...): ...

modified_query_fps: Incomplete

def compute_fraggle_similarity_for_subs(inMol, qMol, qSmi, qSubs, tverskyThresh: float = ...): ...
def GetFraggleSimilarity(queryMol, refMol, tverskyThresh: float = ...): ...
