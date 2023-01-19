from _typeshed import Incomplete

def setDescriptorVersion(version: str = ...): ...

class VectorDescriptorNamespace(dict):
    def __init__(self, **kwargs) -> None: ...

class VectorDescriptorWrapper:
    func: Incomplete
    names: Incomplete
    func_key: Incomplete
    namespace: Incomplete
    def __init__(self, func, names, version, namespace) -> None: ...
    def call_desc(self, mol, index): ...
