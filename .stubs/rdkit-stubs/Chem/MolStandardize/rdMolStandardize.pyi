from typing import Any, ClassVar

from typing import overload
import Boost.Python

class AllowedAtomsValidation(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def validate(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ChargeCorrection(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    Charge: Any
    Name: Any
    Smarts: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class CleanupParameters(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    acidbaseFile: Any
    doCanonical: Any
    fragmentFile: Any
    largestFragmentChooserCountHeavyAtomsOnly: Any
    largestFragmentChooserUseAtomCount: Any
    maxRestarts: Any
    maxTautomers: Any
    maxTransforms: Any
    normalizationsFile: Any
    preferOrganic: Any
    tautomerReassignStereo: Any
    tautomerRemoveBondStereo: Any
    tautomerRemoveIsotopicHs: Any
    tautomerRemoveSp3Stereo: Any
    tautomerTransformsFile: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class DisallowedAtomsValidation(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def validate(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragmentRemover(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def remove(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragmentValidation(MolVSValidations):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def run(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class IsotopeValidation(MolVSValidations):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def run(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class LargestFragmentChooser(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def choose(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MetalDisconnector(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def Disconnect(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetMetalNof(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetMetalNon(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def MetalNof(self) -> Any: ...
    @property
    def MetalNon(self) -> Any: ...

class MolVSValidation(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def validate(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolVSValidations(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def run(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class NeutralValidation(MolVSValidations):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def run(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class NoAtomValidation(MolVSValidations):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def run(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class Normalizer(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def normalize(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RDKitValidation(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def validate(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class Reionizer(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def reionize(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SmilesTautomerMap(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def items(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def keys(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def values(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class Tautomer(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def kekulized(self) -> Any: ...
    @property
    def tautomer(self) -> Any: ...

class TautomerEnumerator(Boost.Python.instance):
    @overload
    @classmethod
    def __init__(cls, classboost) -> Any: ...
    @overload
    @classmethod
    def __init__(cls, classboost, structRDKit) -> Any: ...
    @overload
    @classmethod
    def __init__(cls, classboost, classRDKit) -> Any: ...
    @classmethod
    def Canonicalize(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def Enumerate(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetCallback(cls, classRDKit) -> Any: ...
    @classmethod
    def GetMaxTautomers(cls, classRDKit) -> Any: ...
    @classmethod
    def GetMaxTransforms(cls, classRDKit) -> Any: ...
    @classmethod
    def GetReassignStereo(cls, classRDKit) -> Any: ...
    @classmethod
    def GetRemoveBondStereo(cls, classRDKit) -> Any: ...
    @classmethod
    def GetRemoveSp3Stereo(cls, classRDKit) -> Any: ...
    @classmethod
    def PickCanonical(cls, classRDKit, classboost) -> Any: ...
    def ScoreTautomer(self, *args, **kwargs) -> Any: ...
    @classmethod
    def SetCallback(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetMaxTautomers(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def SetMaxTransforms(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def SetReassignStereo(cls, classRDKit, bool) -> Any: ...
    @classmethod
    def SetRemoveBondStereo(cls, classRDKit, bool) -> Any: ...
    @classmethod
    def SetRemoveSp3Stereo(cls, classRDKit, bool) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def tautomerScoreVersion(self) -> Any: ...

class TautomerEnumeratorCallback(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TautomerEnumeratorResult(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __call__(cls, class `anonymousnamespace') -> Any: ...
    @classmethod
    def __getitem__(cls, class `anonymousnamespace', int) -> Any: ...
    @classmethod
    def __iter__(cls, structboost) -> Any: ...
    @classmethod
    def __len__(cls, class `anonymousnamespace') -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def modifiedAtoms(self) -> Any: ...
    @property
    def modifiedBonds(self) -> Any: ...
    @property
    def smiles(self) -> Any: ...
    @property
    def smilesTautomerMap(self) -> Any: ...
    @property
    def status(self) -> Any: ...
    @property
    def tautomers(self) -> Any: ...

class TautomerEnumeratorStatus(Boost.Python.enum):
    Canceled: ClassVar[TautomerEnumeratorStatus] = ...
    Completed: ClassVar[TautomerEnumeratorStatus] = ...
    MaxTautomersReached: ClassVar[TautomerEnumeratorStatus] = ...
    MaxTransformsReached: ClassVar[TautomerEnumeratorStatus] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class Uncharger(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def uncharge(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class map_indexing_suite_SmilesTautomerMap_entry(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def data(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def key(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def CHARGE_CORRECTIONS() -> Any: ...
def CanonicalTautomer(classRDKit) -> Any: ...
def ChargeParent(classRDKit) -> Any: ...
def Cleanup(classRDKit) -> Any: ...
def FragmentParent(classRDKit) -> Any: ...
def FragmentRemoverFromData(*args, **kwargs) -> Any: ...
def GetV1TautomerEnumerator() -> Any: ...
def IsotopeParent(classRDKit) -> Any: ...
def Normalize(classRDKit) -> Any: ...
def NormalizerFromData(*args, **kwargs) -> Any: ...
def NormalizerFromParams(structRDKit) -> Any: ...
def Reionize(classRDKit) -> Any: ...
def ReionizerFromData(*args, **kwargs) -> Any: ...
def RemoveFragments(classRDKit) -> Any: ...
def StandardizeSmiles(*args, **kwargs) -> Any: ...
def StereoParent(classRDKit) -> Any: ...
def SuperParent(classRDKit) -> Any: ...
def TautomerParent(classRDKit) -> Any: ...
def UpdateParamsFromJSON(*args, **kwargs) -> Any: ...
def ValidateSmiles(*args, **kwargs) -> Any: ...
