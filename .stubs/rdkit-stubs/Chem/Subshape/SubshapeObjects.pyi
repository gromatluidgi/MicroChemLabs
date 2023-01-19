from _typeshed import Incomplete

class SkeletonPoint:
    location: Incomplete
    shapeMoments: Incomplete
    shapeDirs: Incomplete
    molFeatures: Incomplete
    featmapFeatures: Incomplete
    fracVol: float
    def __init__(self, *args, **kwargs) -> None: ...

class ShapeWithSkeleton:
    grid: Incomplete
    skelPts: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...

class SubshapeShape:
    shapes: Incomplete
    featMap: Incomplete
    keyFeat: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...

def DisplaySubshapeSkeleton(viewer, shape, name, color=..., colorByOrder: bool = ...) -> None: ...
def DisplaySubshape(viewer, shape, name, showSkelPts: bool = ..., color=...) -> None: ...
