from typing import Any, ClassVar

import Boost.Python

class FreeChemicalFeature(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetFamily(cls, classChemicalFeatures) -> Any: ...
    @classmethod
    def GetId(cls, classChemicalFeatures) -> Any: ...
    @classmethod
    def GetPos(cls, classChemicalFeatures) -> Any: ...
    @classmethod
    def GetType(cls, classChemicalFeatures) -> Any: ...
    @classmethod
    def SetFamily(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetId(cls, classChemicalFeatures, int) -> Any: ...
    @classmethod
    def SetPos(cls, classChemicalFeatures, classRDGeom) -> Any: ...
    @classmethod
    def SetType(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __getinitargs__(cls, classChemicalFeatures) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
