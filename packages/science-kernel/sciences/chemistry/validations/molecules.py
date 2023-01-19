"""
Molecule Validation Engine

References:

An open source chemical structure curation pipeline using RDKit
https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00456-1#Sec2
"""

from abc import ABC

__RULES = []


class MolValidationRule(ABC):
    """MolValidationRule"""


class MolValidationEngine:
    """MolValidationEngine"""
