from typing import Any, ClassVar

from typing import overload
import Boost.Python

class CachedMolHolder(MolHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddBinary(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class CachedSmilesMolHolder(MolHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddSmiles(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class CachedTrustedSmilesMolHolder(MolHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddSmiles(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FPHolderBase(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddFingerprint(cls, classRDKit, classExplicitBitVect) -> Any: ...
    @classmethod
    def AddMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetFingerprint(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def MakeFingerprint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def PassesFilter(cls, classRDKit, unsignedint, classExplicitBitVect) -> Any: ...
    @classmethod
    def __len__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class KeyFromPropHolder(KeyHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetPropName(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class KeyHolderBase(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddKey(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def AddMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetKey(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def GetKeys(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __len__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolHolder(MolHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolHolderBase(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMol(cls, classRDKit, unsignedint) -> Any: ...
    @overload
    @classmethod
    def __len__(cls, classRDKit) -> Any: ...
    @overload
    @classmethod
    def __len__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class PatternHolder(FPHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SubstructLibrary(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def CountMatches(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetFpHolder(cls, classRDKit) -> Any: ...
    @classmethod
    def GetKeyHolder(cls, classRDKit) -> Any: ...
    @classmethod
    def GetMatches(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMol(cls, classRDKit, unsignedint) -> Any: ...
    @classmethod
    def GetMolHolder(cls, classRDKit) -> Any: ...
    @classmethod
    def GetSearchOrder(cls, classRDKit) -> Any: ...
    @classmethod
    def HasMatch(cls, *args, **kwargs) -> Any: ...
    @overload
    @classmethod
    def InitFromStream(cls, stream) -> Any: ...
    @overload
    @classmethod
    def InitFromStream(cls, f) -> Any: ...
    @overload
    @classmethod
    def InitFromStream(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def Serialize(cls, classRDKit) -> Any: ...
    @classmethod
    def SetSearchOrder(cls, classRDKit, classboost) -> Any: ...
    @overload
    @classmethod
    def ToStream(cls, stream) -> Any: ...
    @overload
    @classmethod
    def ToStream(cls, stream) -> Any: ...
    @overload
    @classmethod
    def ToStream(cls, classRDKit, classboost) -> Any: ...
    @classmethod
    def __getinitargs__(cls, classRDKit) -> Any: ...
    @classmethod
    def __len__(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TautomerPatternHolder(FPHolderBase):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

@overload
def AddPatterns(classRDKit) -> Any: ...
@overload
def AddPatterns(classRDKit, classboost) -> Any: ...
def SubstructLibraryCanSerialize(*args, **kwargs) -> Any: ...
