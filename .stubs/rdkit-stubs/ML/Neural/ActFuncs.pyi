from _typeshed import Incomplete

class ActFunc:
    def __call__(self, x): ...

class Sigmoid(ActFunc):
    def Eval(self, x): ...
    def Deriv(self, x): ...
    def DerivFromVal(self, val): ...
    beta: Incomplete
    def __init__(self, beta: float = ...) -> None: ...

class TanH(ActFunc):
    def Eval(self, x): ...
    def Deriv(self, x): ...
    def DerivFromVal(self, val): ...
    beta: Incomplete
    def __init__(self, beta: float = ...) -> None: ...
