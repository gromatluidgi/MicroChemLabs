from typing import Any, ClassVar

import Boost.Python

class FragCatGenerator(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddFragsFromMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragCatParams(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetFuncGroup(cls, classRDKit, int) -> Any: ...
    @classmethod
    def GetLowerFragLength(cls, classRDKit) -> Any: ...
    @classmethod
    def GetNumFuncGroups(cls, classRDKit) -> Any: ...
    @classmethod
    def GetTolerance(cls, classRDKit) -> Any: ...
    @classmethod
    def GetTypeString(cls, classRDKit) -> Any: ...
    @classmethod
    def GetUpperFragLength(cls, classRDKit) -> Any: ...
    @classmethod
    def Serialize(cls, classRDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragCatalog(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetBitDescription(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetBitDiscrims(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetBitEntryId(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetBitFuncGroupIds(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetBitOrder(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetCatalogParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetEntryBitId(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetEntryDescription(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetEntryDownIds(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetEntryFuncGroupIds(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetEntryOrder(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetFPLength(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetNumEntries(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def Serialize(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __getinitargs__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragFPGenerator(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetFPForMol(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
