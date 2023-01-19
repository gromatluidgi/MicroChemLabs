from typing import Any, ClassVar

from typing import overload
import Boost.Python
ADJUST_IGNOREALL: AdjustQueryWhichFlags
ADJUST_IGNORECHAINS: AdjustQueryWhichFlags
ADJUST_IGNOREDUMMIES: AdjustQueryWhichFlags
ADJUST_IGNORENONDUMMIES: AdjustQueryWhichFlags
ADJUST_IGNORENONE: AdjustQueryWhichFlags
ADJUST_IGNORERINGS: AdjustQueryWhichFlags
AROMATICITY_CUSTOM: AromaticityModel
AROMATICITY_DEFAULT: AromaticityModel
AROMATICITY_MDL: AromaticityModel
AROMATICITY_RDKIT: AromaticityModel
AROMATICITY_SIMPLE: AromaticityModel
LayeredFingerprint_substructLayers: int
SANITIZE_ADJUSTHS: SanitizeFlags
SANITIZE_ALL: SanitizeFlags
SANITIZE_CLEANUP: SanitizeFlags
SANITIZE_CLEANUPCHIRALITY: SanitizeFlags
SANITIZE_FINDRADICALS: SanitizeFlags
SANITIZE_KEKULIZE: SanitizeFlags
SANITIZE_NONE: SanitizeFlags
SANITIZE_PROPERTIES: SanitizeFlags
SANITIZE_SETAROMATICITY: SanitizeFlags
SANITIZE_SETCONJUGATION: SanitizeFlags
SANITIZE_SETHYBRIDIZATION: SanitizeFlags
SANITIZE_SYMMRINGS: SanitizeFlags
_LayeredFingerprint_version: str
_PatternFingerprint_version: str
_RDKFingerprint_version: str

class AdjustQueryParameters(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    adjustConjugatedFiveRings: Any
    adjustDegree: Any
    adjustDegreeFlags: Any
    adjustHeavyDegree: Any
    adjustHeavyDegreeFlags: Any
    adjustRingChain: Any
    adjustRingChainFlags: Any
    adjustRingCount: Any
    adjustRingCountFlags: Any
    adjustSingleBondsBetweenAromaticAtoms: Any
    adjustSingleBondsToDegreeOneNeighbors: Any
    aromatizeIfPossible: Any
    makeAtomsGeneric: Any
    makeAtomsGenericFlags: Any
    makeBondsGeneric: Any
    makeBondsGenericFlags: Any
    makeDummiesQueries: Any
    setMDLFiveRingAromaticity: Any
    useStereoCareForBonds: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def NoAdjustments(self, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AdjustQueryWhichFlags(Boost.Python.enum):
    ADJUST_IGNOREALL: ClassVar[AdjustQueryWhichFlags] = ...
    ADJUST_IGNORECHAINS: ClassVar[AdjustQueryWhichFlags] = ...
    ADJUST_IGNOREDUMMIES: ClassVar[AdjustQueryWhichFlags] = ...
    ADJUST_IGNORENONDUMMIES: ClassVar[AdjustQueryWhichFlags] = ...
    ADJUST_IGNORENONE: ClassVar[AdjustQueryWhichFlags] = ...
    ADJUST_IGNORERINGS: ClassVar[AdjustQueryWhichFlags] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class AromaticityModel(Boost.Python.enum):
    AROMATICITY_CUSTOM: ClassVar[AromaticityModel] = ...
    AROMATICITY_DEFAULT: ClassVar[AromaticityModel] = ...
    AROMATICITY_MDL: ClassVar[AromaticityModel] = ...
    AROMATICITY_RDKIT: ClassVar[AromaticityModel] = ...
    AROMATICITY_SIMPLE: ClassVar[AromaticityModel] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class MolzipLabel(Boost.Python.enum):
    AtomMapNumber: ClassVar[MolzipLabel] = ...
    AtomType: ClassVar[MolzipLabel] = ...
    FragmentOnBonds: ClassVar[MolzipLabel] = ...
    Isotope: ClassVar[MolzipLabel] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class MolzipParams(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    enforceValenceRules: Any
    label: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def setAtomSymbols(cls, structRDKit, classboost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RemoveHsParameters(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    removeAndTrackIsotopes: Any
    removeDefiningBondStereo: Any
    removeDegreeZero: Any
    removeDummyNeighbors: Any
    removeHigherDegrees: Any
    removeHydrides: Any
    removeInSGroups: Any
    removeIsotopes: Any
    removeMapped: Any
    removeNonimplicit: Any
    removeNontetrahedralNeighbors: Any
    removeOnlyHNeighbors: Any
    removeWithQuery: Any
    removeWithWedgedBond: Any
    showWarnings: Any
    updateExplicitCount: Any
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUSTHS: ClassVar[SanitizeFlags] = ...
    SANITIZE_ALL: ClassVar[SanitizeFlags] = ...
    SANITIZE_CLEANUP: ClassVar[SanitizeFlags] = ...
    SANITIZE_CLEANUPCHIRALITY: ClassVar[SanitizeFlags] = ...
    SANITIZE_FINDRADICALS: ClassVar[SanitizeFlags] = ...
    SANITIZE_KEKULIZE: ClassVar[SanitizeFlags] = ...
    SANITIZE_NONE: ClassVar[SanitizeFlags] = ...
    SANITIZE_PROPERTIES: ClassVar[SanitizeFlags] = ...
    SANITIZE_SETAROMATICITY: ClassVar[SanitizeFlags] = ...
    SANITIZE_SETCONJUGATION: ClassVar[SanitizeFlags] = ...
    SANITIZE_SETHYBRIDIZATION: ClassVar[SanitizeFlags] = ...
    SANITIZE_SYMMRINGS: ClassVar[SanitizeFlags] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class StereoBondThresholds(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def CHIRAL_ATOM(self) -> Any: ...
    @property
    def DBL_BOND_NO_STEREO(self) -> Any: ...
    @property
    def DBL_BOND_SPECIFIED_STEREO(self) -> Any: ...
    @property
    def DIRECTION_SET(self) -> Any: ...

class _vectstruct RDKit::Chirality::StereoInfo(Boost.Python.instance):
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

def AddHs(classRDKit) -> Any: ...
def AddRecursiveQuery(*args, **kwargs) -> Any: ...
def AddWavyBondsForStereoAny(classRDKit) -> Any: ...
def AdjustQueryProperties(classRDKit) -> Any: ...
def AssignAtomChiralTagsFromMolParity(classRDKit) -> Any: ...
def AssignAtomChiralTagsFromStructure(classRDKit) -> Any: ...
def AssignChiralTypesFromBondDirs(classRDKit) -> Any: ...
def AssignRadicals(classRDKit) -> Any: ...
def AssignStereochemistry(*args, **kwargs) -> Any: ...
def AssignStereochemistryFrom3D(*args, **kwargs) -> Any: ...
def Cleanup(classRDKit) -> Any: ...
def CombineMols(*args, **kwargs) -> Any: ...
def ConvertGenericQueriesToSubstanceGroups(classRDKit) -> Any: ...
def DeleteSubstructs(*args, **kwargs) -> Any: ...
def DetectBondStereoChemistry(*args, **kwargs) -> Any: ...
def DetectBondStereochemistry(classRDKit) -> Any: ...
def DetectChemistryProblems(classRDKit) -> Any: ...
def FastFindRings(classRDKit) -> Any: ...
def FindAllPathsOfLengthN(*args, **kwargs) -> Any: ...
def FindAllSubgraphsOfLengthMToN(*args, **kwargs) -> Any: ...
def FindAllSubgraphsOfLengthN(*args, **kwargs) -> Any: ...
def FindAtomEnvironmentOfRadiusN(*args, **kwargs) -> Any: ...
def FindPotentialStereo(classRDKit) -> Any: ...
def FindPotentialStereoBonds(*args, **kwargs) -> Any: ...
def FindRingFamilies(classRDKit) -> Any: ...
def FindUniqueSubgraphsOfLengthN(*args, **kwargs) -> Any: ...
def FragmentOnBRICSBonds(classRDKit) -> Any: ...
def FragmentOnBonds(classRDKit, classboost) -> Any: ...
def FragmentOnSomeBonds(classRDKit, classboost) -> Any: ...
def Get3DDistanceMatrix(classRDKit) -> Any: ...
def GetAdjacencyMatrix(classRDKit) -> Any: ...
def GetAllowNontetrahedralChirality() -> Any: ...
def GetDistanceMatrix(classRDKit) -> Any: ...
def GetFormalCharge(classRDKit) -> Any: ...
@overload
def GetMolFrags() -> Any: ...
@overload
def GetMolFrags(classRDKit) -> Any: ...
def GetMostSubstitutedCoreMatch(*args, **kwargs) -> Any: ...
def GetSSSR(classRDKit) -> Any: ...
def GetShortestPath(*args, **kwargs) -> Any: ...
def GetSymmSSSR(classRDKit) -> Any: ...
def GetUseLegacyStereoPerception() -> Any: ...
def Kekulize(classRDKit) -> Any: ...
def KekulizeIfPossible(classRDKit) -> Any: ...
def LayeredFingerprint(classRDKit) -> Any: ...
def MergeQueryHs(classRDKit) -> Any: ...
def MolAddRecursiveQueries(*args, **kwargs) -> Any: ...
def MurckoDecompose(classRDKit) -> Any: ...
def ParseMolQueryDefFile(classboost) -> Any: ...
def PathToSubmol(classRDKit, classboost) -> Any: ...
@overload
def PatternFingerprint(classRDKit) -> Any: ...
@overload
def PatternFingerprint(classRDKit) -> Any: ...
def RDKFingerprint(*args, **kwargs) -> Any: ...
def ReapplyMolBlockWedging(classRDKit) -> Any: ...
def RemoveAllHs(classRDKit) -> Any: ...
@overload
def RemoveHs(classRDKit) -> Any: ...
@overload
def RemoveHs(classRDKit, structRDKit) -> Any: ...
def RemoveStereochemistry(classRDKit) -> Any: ...
def RenumberAtoms(classRDKit, classboost) -> Any: ...
def ReplaceCore(mol, core, match) -> Any: ...
def ReplaceSidechains(*args, **kwargs) -> Any: ...
def ReplaceSubstructs(*args, **kwargs) -> Any: ...
def SanitizeMol(*args, **kwargs) -> Any: ...
def SetAllowNontetrahedralChirality(bool) -> Any: ...
def SetAromaticity(classRDKit) -> Any: ...
def SetBondStereoFromDirections(classRDKit) -> Any: ...
def SetConjugation(classRDKit) -> Any: ...
def SetDoubleBondNeighborDirections(classRDKit) -> Any: ...
def SetGenericQueriesFromProperties(classRDKit) -> Any: ...
def SetHybridization(classRDKit) -> Any: ...
def SetTerminalAtomCoords(*args, **kwargs) -> Any: ...
def SetUseLegacyStereoPerception(bool) -> Any: ...
def SortMatchesByDegreeOfCoreSubstitution(*args, **kwargs) -> Any: ...
def SplitMolByPDBChainId(classRDKit) -> Any: ...
def SplitMolByPDBResidues(classRDKit) -> Any: ...
def UnfoldedRDKFingerprintCountBased(classRDKit) -> Any: ...
def WedgeBond(*args, **kwargs) -> Any: ...
def WedgeMolBonds(*args, **kwargs) -> Any: ...
def molzip(classRDKit) -> Any: ...