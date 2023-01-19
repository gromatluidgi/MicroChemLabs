from _typeshed import Incomplete
from rdkit import Geometry as Geometry, RDConfig as RDConfig

logger: Incomplete
BEGIN: int
END: int
TRIANGLE_FAN: int
COLOR: int
VERTEX: int
NORMAL: int
SPHERE: int
CYLINDER: int
ALPHA: int

def ShowArrow(viewer, tail, head, radius, color, label, transparency: int = ..., includeArrowhead: bool = ...) -> None: ...
def ShowMolFeats(mol, factory, viewer, radius: float = ..., confId: int = ..., showOnly: bool = ..., name: str = ..., transparency: float = ..., colors: Incomplete | None = ..., excludeTypes=..., useFeatDirs: bool = ..., featLabel: Incomplete | None = ..., dirLabel: Incomplete | None = ..., includeArrowheads: bool = ..., writeFeats: bool = ..., showMol: bool = ..., featMapFile: bool = ...) -> None: ...

parser: Incomplete
