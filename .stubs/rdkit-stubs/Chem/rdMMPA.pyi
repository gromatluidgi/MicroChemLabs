from typing import Any

from typing import overload

@overload
def FragmentMol(classRDKit) -> Any: ...
@overload
def FragmentMol(classRDKit, classboost) -> Any: ...
