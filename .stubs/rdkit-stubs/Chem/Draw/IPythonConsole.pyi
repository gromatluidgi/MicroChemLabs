from . import InteractiveRenderer as InteractiveRenderer
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Draw as Draw, rdChemReactions as rdChemReactions, rdchem as rdchem
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D

molSize: Incomplete
highlightSubstructs: bool
kekulizeStructures: bool
highlightByReactant: bool
ipython_useSVG: bool
ipython_showProperties: bool
ipython_maxProperties: int
ipython_3d: bool
molSize_3d: Incomplete
drawing_type_3d: str
bgcolor_3d: str
drawOptions: Incomplete

def addMolToView(mol, view, confId: int = ..., drawAs: Incomplete | None = ...) -> None: ...
def drawMol3D(m, view: Incomplete | None = ..., confId: int = ..., drawAs: Incomplete | None = ..., bgColor: Incomplete | None = ..., size: Incomplete | None = ...): ...
def display_pil_image(img): ...
def ShowMols(mols, maxMols: int = ..., **kwargs): ...
def DrawMorganBit(mol, bitId, bitInfo, drawOptions=..., **kwargs): ...
def DrawMorganBits(*args, drawOptions=..., **kwargs): ...
def DrawRDKitBit(mol, bitId, bitInfo, drawOptions=..., **kwargs): ...
def DrawRDKitBits(*args, drawOptions=..., **kwargs): ...
def EnableSubstructMatchRendering() -> None: ...
def InstallIPythonRenderer(): ...
def DisableSubstructMatchRendering() -> None: ...
def UninstallIPythonRenderer() -> None: ...