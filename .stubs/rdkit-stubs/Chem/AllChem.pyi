from rdkit.Chem import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdDistGeom import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdReducedGraphs import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Chem.rdqueries import *
from rdkit.Chem.rdMolEnumerator import *
from rdkit.Chem.rdSLNParse import *
from _typeshed import Incomplete
from collections.abc import Generator
from rdkit import DataStructs as DataStructs, ForceField as ForceField, RDConfig as RDConfig, rdBase as rdBase
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers as EnumerateStereoisomers, StereoEnumerationOptions as StereoEnumerationOptions
from rdkit.Geometry import rdGeometry as rdGeometry
from rdkit.RDLogger import logger as logger

def TransformMol(mol, tform, confId: int = ..., keepConfs: bool = ...) -> None: ...
def ComputeMolShape(mol, confId: int = ..., boxDim=..., spacing: float = ..., **kwargs): ...
def ComputeMolVolume(mol, confId: int = ..., gridSpacing: float = ..., boxMargin: float = ...): ...
def GetConformerRMS(mol, confId1, confId2, atomIds: Incomplete | None = ..., prealigned: bool = ...): ...
def GetConformerRMSMatrix(mol, atomIds: Incomplete | None = ..., prealigned: bool = ...): ...
def EnumerateLibraryFromReaction(reaction, sidechainSets, returnReactants: bool = ...) -> Generator[Incomplete, None, None]: ...
def ConstrainedEmbed(mol, core, useTethers: bool = ..., coreConfId: int = ..., randomseed: int = ..., getForceField=..., **kwargs): ...
def AssignBondOrdersFromTemplate(refmol, mol): ...
