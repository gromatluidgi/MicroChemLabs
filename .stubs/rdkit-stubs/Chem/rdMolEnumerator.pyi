from typing import Any, ClassVar

from typing import overload
import Boost.Python

class EnumeratorType(Boost.Python.enum):
    LinkNode: ClassVar[EnumeratorType] = ...
    PositionVariation: ClassVar[EnumeratorType] = ...
    RepeatUnit: ClassVar[EnumeratorType] = ...
    names: ClassVar[dict] = ...
    values: ClassVar[dict] = ...
    __slots__: ClassVar[tuple] = ...

class MolEnumeratorParams(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...
    doRandom: Any
    maxToEnumerate: Any
    randomSeed: Any
    sanitize: Any
    @classmethod
    def __init__(cls, classboost, enum `anonymousnamespace') -> Any: ...
    @classmethod
    def SetEnumerationOperator(cls, structRDKit, enum `anonymousnamespace') -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

@overload
def Enumerate(classRDKit) -> Any: ...
@overload
def Enumerate(classRDKit, structRDKit) -> Any: ...
