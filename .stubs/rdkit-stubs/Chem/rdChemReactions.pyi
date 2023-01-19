from typing import Any, ClassVar

from typing import overload
import Boost.Python
SANITIZE_ADJUST_REACTANTS: SanitizeFlags
SANITIZE_ALL: SanitizeFlags
SANITIZE_ATOM_MAPS: SanitizeFlags
SANITIZE_MERGEHS: SanitizeFlags
SANITIZE_NONE: SanitizeFlags
SANITIZE_RGROUP_NAMES: SanitizeFlags

class CartesianProductStrategy(EnumerationStrategyBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ChemicalReaction(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddAgentTemplate(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def AddProductTemplate(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def AddReactantTemplate(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def AddRecursiveQueriesToReaction(cls, classRDKit) -> Any: ...
    @classmethod
    def ClearComputedProps(cls, classRDKit) -> Any: ...
    @classmethod
    def ClearProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetAgentTemplate(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def GetAgents(cls, classRDKit) -> Any: ...
    @classmethod
    def GetBoolProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetDoubleProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetIntProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetNumAgentTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def GetNumProductTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def GetNumReactantTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def GetProductTemplate(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def GetProducts(cls, classRDKit) -> Any: ...
    @classmethod
    def GetProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetPropNames(cls, classRDKit) -> Any: ...
    @classmethod
    def GetPropsAsDict(cls, classRDKit) -> Any: ...
    @classmethod
    def GetReactantTemplate(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def GetReactants(cls, classRDKit) -> Any: ...
    @classmethod
    def GetReactingAtoms(cls, classRDKit) -> Any: ...
    @classmethod
    def GetUnsignedProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def HasProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def Initialize(cls, classRDKit) -> Any: ...
    @classmethod
    def IsInitialized(cls, classRDKit) -> Any: ...
    @classmethod
    def IsMoleculeAgent(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def IsMoleculeProduct(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def IsMoleculeReactant(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def RemoveAgentTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def RemoveUnmappedProductTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def RemoveUnmappedReactantTemplates(cls, classRDKit) -> Any: ...
    @classmethod
    def RunReactant(cls, classRDKit, classboost, unsignedint) -> Any: ...
    @classmethod
    def RunReactantInPlace(cls, *args, **kwargs) -> Any: ...
    @overload
    @classmethod
    def RunReactants(cls, classRDKit, classboost) -> Any: ...
    @overload
    @classmethod
    def RunReactants(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def SetBoolProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetDoubleProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetIntProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetProp(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetUnsignedProp(cls, *args, **kwargs) -> Any: ...
    @overload
    @classmethod
    def ToBinary(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def ToBinary(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def Validate(cls, classRDKit) -> Any: ...
    @classmethod
    def _getImplicitPropertiesFlag(cls, classRDKit) -> Any: ...
    @classmethod
    def _setImplicitPropertiesFlag(cls, classRDKit, bool) -> Any: ...
    @classmethod
    def __getinitargs__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerateLibrary(EnumerateLibraryBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetReagents(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerateLibraryBase(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetEnumerator(cls, classRDKit) -> Any: ...
    @classmethod
    def GetPosition(cls, classRDKit) -> Any: ...
    @classmethod
    def GetReaction(cls, classRDKit) -> Any: ...
    @classmethod
    def GetState(cls, classRDKit) -> Any: ...
    @classmethod
    def InitFromString(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def ResetState(cls, classRDKit) -> Any: ...
    @classmethod
    def Serialize(cls, classRDKit) -> Any: ...
    @classmethod
    def SetState(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def next(cls, classRDKit) -> Any: ...
    @classmethod
    def nextSmiles(cls, classRDKit) -> Any: ...
    @classmethod
    def __bool__(cls, classRDKit) -> Any: ...
    @classmethod
    def __iter__(cls, classboost) -> Any: ...
    @classmethod
    def __next__(cls, classRDKit) -> Any: ...
    @classmethod
    def __nonzero__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerationParams(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    reagentMaxMatchCount: Any
    sanePartialProducts: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EnumerationStrategyBase(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetNumPermutations(cls, classRDKit) -> Any: ...
    @classmethod
    def GetPosition(cls, classRDKit) -> Any: ...
    @classmethod
    def Initialize(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def Skip(cls, classRDKit, unsigned__int64) -> Any: ...
    @classmethod
    def Type(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def next(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def next(cls, classRDKit) -> Any: ...
    @classmethod
    def __bool__(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def __next__(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def __next__(cls, classRDKit) -> Any: ...
    @classmethod
    def __nonzero__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class EvenSamplePairsStrategy(EnumerationStrategyBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def Stats(cls, classRDKit) -> Any: ...
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FingerprintType(Boost.Python.enum):
    AtomPairFP: ClassVar[FingerprintType] = ...
    MorganFP: ClassVar[FingerprintType] = ...
    PatternFP: ClassVar[FingerprintType] = ...
    RDKitFP: ClassVar[FingerprintType] = ...
    TopologicalTorsion: ClassVar[FingerprintType] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class MOL_SPTR_VECT(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def append(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, structboost, classstd) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class RandomSampleAllBBsStrategy(EnumerationStrategyBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RandomSampleStrategy(EnumerationStrategyBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __copy__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class ReactionFingerprintParams(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    agentWeight: Any
    bitRatioAgents: Any
    fpSize: Any
    fpType: Any
    includeAgents: Any
    nonAgentWeight: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUST_REACTANTS: ClassVar[SanitizeFlags] = ...
    SANITIZE_ALL: ClassVar[SanitizeFlags] = ...
    SANITIZE_ATOM_MAPS: ClassVar[SanitizeFlags] = ...
    SANITIZE_MERGEHS: ClassVar[SanitizeFlags] = ...
    SANITIZE_NONE: ClassVar[SanitizeFlags] = ...
    SANITIZE_RGROUP_NAMES: ClassVar[SanitizeFlags] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class VectMolVect(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def append(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(cls, *args, **kwargs) -> Any: ...
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

class VectSizeT(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def append(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, structboost, classstd) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

class VectorOfStringVectors(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def append(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(cls, *args, **kwargs) -> Any: ...
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

def Compute2DCoordsForReaction(classRDKit) -> Any: ...
def CreateDifferenceFingerprintForReaction(classRDKit) -> Any: ...
def CreateStructuralFingerprintForReaction(classRDKit) -> Any: ...
def EnumerateLibraryCanSerialize(*args, **kwargs) -> Any: ...
def GetChemDrawRxnAdjustParams() -> Any: ...
def GetDefaultAdjustParams() -> Any: ...
def HasAgentTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasProductTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasReactantTemplateSubstructMatch(*args, **kwargs) -> Any: ...
def HasReactionAtomMapping(classRDKit) -> Any: ...
def HasReactionSubstructMatch(*args, **kwargs) -> Any: ...
def IsReactionTemplateMoleculeAgent(classRDKit, double) -> Any: ...
def MatchOnlyAtRgroupsAdjustParams() -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(rxn) -> Any: ...
@overload
def PreprocessReaction(classRDKit) -> Any: ...
def ReactionFromMolecule(classRDKit) -> Any: ...
def ReactionFromPNGFile(*args, **kwargs) -> Any: ...
def ReactionFromPNGString(*args, **kwargs) -> Any: ...
def ReactionFromRxnBlock(*args, **kwargs) -> Any: ...
def ReactionFromRxnFile(*args, **kwargs) -> Any: ...
def ReactionFromSmarts(*args, **kwargs) -> Any: ...
def ReactionMetadataToPNGFile(classRDKit, classboost) -> Any: ...
def ReactionMetadataToPNGString(classRDKit, classboost) -> Any: ...
def ReactionToMolecule(classRDKit) -> Any: ...
def ReactionToRxnBlock(classRDKit) -> Any: ...
def ReactionToSmarts(classRDKit) -> Any: ...
def ReactionToSmiles(classRDKit) -> Any: ...
def ReactionToV3KRxnBlock(classRDKit) -> Any: ...
def ReduceProductToSideChains(classboost) -> Any: ...
def RemoveMappingNumbersFromReactions(classRDKit) -> Any: ...
def SanitizeRxn(*args, **kwargs) -> Any: ...
def UpdateProductsStereochemistry(classRDKit) -> Any: ...