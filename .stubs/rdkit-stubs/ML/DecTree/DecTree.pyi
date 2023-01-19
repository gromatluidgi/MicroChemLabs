from _typeshed import Incomplete
from rdkit.ML.DecTree import Tree as Tree

class DecTreeNode(Tree.TreeNode):
    examples: Incomplete
    badExamples: Incomplete
    trainingExamples: Incomplete
    testExamples: Incomplete
    def __init__(self, *args, **kwargs) -> None: ...
    def ClassifyExample(self, example, appendExamples: int = ...): ...
    def AddChild(self, name, label: Incomplete | None = ..., data: Incomplete | None = ..., isTerminal: int = ...): ...
    def GetExamples(self): ...
    def SetExamples(self, examples) -> None: ...
    def GetBadExamples(self): ...
    def SetBadExamples(self, examples) -> None: ...
    def GetTrainingExamples(self): ...
    def SetTrainingExamples(self, examples) -> None: ...
    def GetTestExamples(self): ...
    def SetTestExamples(self, examples) -> None: ...
    def ClearExamples(self) -> None: ...