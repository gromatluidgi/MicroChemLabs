from typing import Any, ClassVar

import Boost.Python
avalonSSSBits: int
avalonSimilarityBits: int

class StruChkFlag(Boost.Python.enum):
    alias_conversion_failed: ClassVar[StruChkFlag] = ...
    atom_check_failed: ClassVar[StruChkFlag] = ...
    atom_clash: ClassVar[StruChkFlag] = ...
    bad_molecule: ClassVar[StruChkFlag] = ...
    dubious_stereo_removed: ClassVar[StruChkFlag] = ...
    either_warning: ClassVar[StruChkFlag] = ...
    fragments_found: ClassVar[StruChkFlag] = ...
    names: ClassVar[dict] = ...
    recharged: ClassVar[StruChkFlag] = ...
    size_check_failed: ClassVar[StruChkFlag] = ...
    stereo_error: ClassVar[StruChkFlag] = ...
    stereo_forced_bad: ClassVar[StruChkFlag] = ...
    stereo_transformed: ClassVar[StruChkFlag] = ...
    template_transformed: ClassVar[StruChkFlag] = ...
    transformed: ClassVar[StruChkFlag] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class StruChkResult(Boost.Python.enum):
    bad_set: ClassVar[StruChkResult] = ...
    names: ClassVar[dict] = ...
    success: ClassVar[StruChkResult] = ...
    transformed_set: ClassVar[StruChkResult] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

def CheckMolecule(classRDKit) -> Any: ...
def CheckMoleculeString(*args, **kwargs) -> Any: ...
def CloseCheckMolFiles() -> Any: ...
def Generate2DCoords(classRDKit) -> Any: ...
def GetAvalonCountFP(classRDKit) -> Any: ...
def GetAvalonFP(classRDKit) -> Any: ...
def GetAvalonFPAsWords(classRDKit) -> Any: ...
def GetCanonSmiles(classRDKit) -> Any: ...
def GetCheckMolLog() -> Any: ...
def InitializeCheckMol(*args, **kwargs) -> Any: ...