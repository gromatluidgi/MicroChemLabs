from . import ClusterUtils as ClusterUtils
from _typeshed import Incomplete
from rdkit.piddle import piddle as piddle
from rdkit.sping import pid as pid

piddle = pid

class VisOpts:
    xOffset: int
    yOffset: int
    lineColor: Incomplete
    hideColor: Incomplete
    terminalColors: Incomplete
    lineWidth: int
    hideWidth: float
    nodeRad: int
    nodeColor: Incomplete
    highlightColor: Incomplete
    highlightRad: int

class ClusterRenderer:
    canvas: Incomplete
    size: Incomplete
    ptColors: Incomplete
    lineWidth: Incomplete
    showIndices: Incomplete
    showNodes: Incomplete
    stopAtCentroids: Incomplete
    logScale: Incomplete
    tooClose: Incomplete
    def __init__(self, canvas, size, ptColors=..., lineWidth: Incomplete | None = ..., showIndices: int = ..., showNodes: int = ..., stopAtCentroids: int = ..., logScale: int = ..., tooClose: int = ...) -> None: ...
    ySpace: Incomplete
    def DrawTree(self, cluster, minHeight: float = ...) -> None: ...

def DrawClusterTree(cluster, canvas, size, ptColors=..., lineWidth: Incomplete | None = ..., showIndices: int = ..., showNodes: int = ..., stopAtCentroids: int = ..., logScale: int = ..., tooClose: int = ...) -> None: ...
def ClusterToPDF(cluster, fileName, size=..., ptColors=..., lineWidth: Incomplete | None = ..., showIndices: int = ..., stopAtCentroids: int = ..., logScale: int = ...): ...
def ClusterToSVG(cluster, fileName, size=..., ptColors=..., lineWidth: Incomplete | None = ..., showIndices: int = ..., stopAtCentroids: int = ..., logScale: int = ...): ...
def ClusterToImg(cluster, fileName, size=..., ptColors=..., lineWidth: Incomplete | None = ..., showIndices: int = ..., stopAtCentroids: int = ..., logScale: int = ...): ...
