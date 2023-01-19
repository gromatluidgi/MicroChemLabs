from typing import Any, ClassVar

from typing import overload
import Boost.Python

class EmbedParameters(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    ETversion: Any
    boundsMatForceScaling: Any
    boxSizeMult: Any
    clearConfs: Any
    embedFragmentsSeparately: Any
    enforceChirality: Any
    forceTransAmides: Any
    ignoreSmoothingFailures: Any
    maxIterations: Any
    numThreads: Any
    numZeroFail: Any
    onlyHeavyAtomsForRMS: Any
    optimizerForceTol: Any
    pruneRmsThresh: used to filter multiple conformations
    randNegEig: Any
    randomSeed: Any
    useBasicKnowledge: Any
    useExpTorsionAnglePrefs: Any
    useMacrocycleTorsions: Any
    useRandomCoords: Any
    useSmallRingTorsions: Any
    useSymmetryForPruning: Any
    verbose: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def SetBoundsMat(cls, structRDKit, classboost) -> Any: ...
    @classmethod
    def SetCPCI(cls, structRDKit, classboost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def ETDG() -> Any: ...
def ETKDG() -> Any: ...
def ETKDGv2() -> Any: ...
def ETKDGv3() -> Any: ...
@overload
def EmbedMolecule(classRDKit) -> Any: ...
@overload
def EmbedMolecule(classRDKit, structRDKit) -> Any: ...
def EmbedMultipleConfs(*args, **kwargs) -> Any: ...
def GetMoleculeBoundsMatrix(classRDKit) -> Any: ...
def KDG() -> Any: ...
def srETKDGv3() -> Any: ...
