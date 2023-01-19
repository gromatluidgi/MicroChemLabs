from _typeshed import Incomplete
from rdkit.Chem.Draw.canvasbase import CanvasBase as CanvasBase
from rdkit.Chem.Draw.rdMolDraw2DQt import rdkitQtVersion as rdkitQtVersion

QPainter_Antialiasing: Incomplete
QPainter_SmoothPixmapTransform: Incomplete
GlobalColor_white: Incomplete
PenStile_SolidLine: Incomplete
PenStile_DashLine: Incomplete

class Canvas(CanvasBase):
    size: Incomplete
    qsize: Incomplete
    pixmap: Incomplete
    painter: Incomplete
    def __init__(self, size) -> None: ...
    def addCanvasLine(self, p1, p2, color=..., color2: Incomplete | None = ..., **kwargs) -> None: ...
    def addCanvasText(self, text, pos, font, color=..., **kwargs): ...
    def addCanvasPolygon(self, ps, color=..., fill: bool = ..., stroke: bool = ..., **kwargs) -> None: ...
    def addCanvasDashedWedge(self, p1, p2, p3, dash=..., color=..., color2: Incomplete | None = ..., **kwargs) -> None: ...
    def flush(self) -> None: ...
