from typing import Any

import Boost.Python

class ForceField(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def AddDistanceConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def AddExtraPoint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def AddFixedPoint(cls, classForceFields, unsignedint) -> Any: ...
    @classmethod
    def CalcEnergy(cls, classForceFields) -> Any: ...
    @classmethod
    def CalcGrad(cls, classForceFields) -> Any: ...
    @classmethod
    def Dimension(cls, classForceFields) -> Any: ...
    @classmethod
    def GetExtraPointPos(cls, classForceFields, unsignedint) -> Any: ...
    @classmethod
    def Initialize(cls, classForceFields) -> Any: ...
    @classmethod
    def MMFFAddAngleConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddDistanceConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddPositionConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def MMFFAddTorsionConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def Minimize(cls, classForceFields) -> Any: ...
    @classmethod
    def MinimizeTrajectory(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def NumPoints(cls, classForceFields) -> Any: ...
    @classmethod
    def Positions(cls, classForceFields) -> Any: ...
    @classmethod
    def UFFAddAngleConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddDistanceConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddPositionConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def UFFAddTorsionConstraint(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MMFFMolProperties(Boost.Python.instance):
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def GetMMFFAngleBendParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFAtomType(cls, classForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFBondStretchParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFFormalCharge(cls, classForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFOopBendParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFPartialCharge(cls, classForceFields, unsignedint) -> Any: ...
    @classmethod
    def GetMMFFStretchBendParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFTorsionParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def GetMMFFVdWParams(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def SetMMFFAngleTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFBondTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFDielectricConstant(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFDielectricModel(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFEleTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFOopTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFStretchBendTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFTorsionTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFVariant(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFVdWTerm(cls, classForceFields) -> Any: ...
    @classmethod
    def SetMMFFVerbosity(cls, classForceFields) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
