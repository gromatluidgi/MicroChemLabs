from qt import *
from _typeshed import Incomplete
from rdkit.qtGui.qtActiveX import MakeActiveXClass as MakeActiveXClass

cdxModule: Incomplete

class ChemdrawPanel(QWidget):
    cdx: Incomplete
    offset: int
    label: Incomplete
    def __init__(self, parent: Incomplete | None = ..., name: str = ..., readOnly: int = ..., size=...) -> None: ...
    def pullData(self, fmt: str = ...): ...
    def setData(self, data, fmt: str = ...): ...
    def resizeEvent(self, evt) -> None: ...
    def __del__(self) -> None: ...
