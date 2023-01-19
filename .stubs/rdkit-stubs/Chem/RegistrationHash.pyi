import enum
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import EnumerateStereoisomers as EnumerateStereoisomers, rdMolHash as rdMolHash
from typing import Iterable, Optional

ATOM_PROP_MAP_NUMBER: str
logger: Incomplete
ENHANCED_STEREO_GROUP_REGEX: Incomplete
ENHANCED_STEREO_GROUP_WEIGHTS: Incomplete
EMPTY_MOL_TAUTOMER_HASH: str

class HashLayer(enum.Enum):
    CANONICAL_SMILES: Incomplete
    ESCAPE: Incomplete
    FORMULA: Incomplete
    NO_STEREO_SMILES: Incomplete
    NO_STEREO_TAUTOMER_HASH: Incomplete
    SGROUP_DATA: Incomplete
    TAUTOMER_HASH: Incomplete

class HashScheme(enum.Enum):
    ALL_LAYERS: Incomplete
    STEREO_INSENSITIVE_LAYERS: Incomplete
    TAUTOMER_INSENSITIVE_LAYERS: Incomplete

def GetMolHash(all_layers, hash_scheme: HashScheme = ...) -> str: ...
def GetMolLayers(original_molecule: Chem.rdchem.Mol, data_field_names: Optional[Iterable] = ..., escape: Optional[str] = ...) -> None: ...
def GetStereoTautomerHash(molecule): ...
def GetCanonicalSmiles(cxsmiles): ...
def GetNoStereoLayers(mol): ...

class EnhancedStereoUpdateMode(enum.Enum):
    ADD_WEIGHTS: Incomplete
    REMOVE_WEIGHTS: Incomplete
