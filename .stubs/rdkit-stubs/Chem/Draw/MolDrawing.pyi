from _typeshed import Incomplete
from rdkit import Chem as Chem

def cmp(t1, t2): ...

periodicTable: Incomplete

class Font:
    face: Incomplete
    size: Incomplete
    weight: Incomplete
    name: Incomplete
    def __init__(self, face: Incomplete | None = ..., size: Incomplete | None = ..., name: Incomplete | None = ..., weight: Incomplete | None = ...) -> None: ...

class DrawingOptions:
    dotsPerAngstrom: int
    useFraction: float
    atomLabelFontFace: str
    atomLabelFontSize: int
    atomLabelMinFontSize: int
    atomLabelDeuteriumTritium: bool
    bondLineWidth: float
    dblBondOffset: float
    dblBondLengthFrac: float
    defaultColor: Incomplete
    selectColor: Incomplete
    bgColor: Incomplete
    colorBonds: bool
    noCarbonSymbols: bool
    includeAtomNumbers: bool
    atomNumberOffset: int
    radicalSymbol: str
    dash: Incomplete
    wedgeDashedBonds: bool
    showUnknownDoubleBonds: bool
    coordScale: float
    elemDict: Incomplete

class MolDrawing:
    canvas: Incomplete
    canvasSize: Incomplete
    drawingOptions: Incomplete
    atomPs: Incomplete
    boundingBoxes: Incomplete
    def __init__(self, canvas: Incomplete | None = ..., drawingOptions: Incomplete | None = ...) -> None: ...
    def transformPoint(self, pos): ...
    molTrans: Incomplete
    currAtomLabelFontSize: Incomplete
    drawingTrans: Incomplete
    def scaleAndCenter(self, mol, conf, coordCenter: bool = ..., canvasSize: Incomplete | None = ..., ignoreHs: bool = ...) -> None: ...
    currDotsPerAngstrom: Incomplete
    activeMol: Incomplete
    bondRings: Incomplete
    def AddMol(self, mol, centerIt: bool = ..., molTrans: Incomplete | None = ..., drawingTrans: Incomplete | None = ..., highlightAtoms=..., confId: int = ..., flagCloseContactsDist: int = ..., highlightMap: Incomplete | None = ..., ignoreHs: bool = ..., highlightBonds=..., **kwargs) -> None: ...
