from _typeshed import Incomplete
from rdkit import Chem as Chem, RDConfig as RDConfig
from rdkit.Avalon import pyAvalonTools as pyAvalonTools
from rdkit.Chem.MolKey import InchiInfo as InchiInfo
from typing import NamedTuple

class MolIdentifierException(Exception): ...
class BadMoleculeException(Exception): ...

MOL_KEY_VERSION: str
ERROR_DICT: Incomplete
INCHI_COMPUTATION_ERROR: Incomplete
RDKIT_CONVERSION_ERROR: Incomplete
INCHI_READWRITE_ERROR: Incomplete
NULL_MOL: Incomplete
BAD_SET: Incomplete
GET_STEREO_RE: Incomplete
NULL_SMILES_RE: Incomplete
PATTERN_NULL_MOL: str
CHIRAL_POS: int
T_NULL_MOL: Incomplete
stereo_code_dict: Incomplete

def initStruchk(configDir: Incomplete | None = ..., logFile: Incomplete | None = ...) -> None: ...
def CheckCTAB(ctab, isSmiles: bool = ...): ...

class InchiResult(NamedTuple):
    error: Incomplete
    inchi: Incomplete
    fixed_ctab: Incomplete

def GetInchiForCTAB(ctab): ...
def ErrorBitsToText(err): ...

class MolKeyResult(NamedTuple):
    mol_key: Incomplete
    error: Incomplete
    inchi: Incomplete
    fixed_ctab: Incomplete
    stereo_code: Incomplete
    stereo_comment: Incomplete

def GetKeyForCTAB(ctab, stereo_info: Incomplete | None = ..., stereo_comment: Incomplete | None = ..., logger: Incomplete | None = ...): ...
