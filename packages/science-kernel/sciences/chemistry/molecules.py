from chembl_structure_pipeline import standardizer
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Mol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")


class MoleculeUtils:
    """Miscellaneous wrapper function for RDKit and other scientific libraries."""

    @staticmethod
    def has_explicit_hygrogen(mol: Mol) -> bool:
        """Indicates if the molecule has explicit hydrogen defined.

        This function can be used to determine if a molecule generated from an SDF
        has a 2D or 3D representation.
        """
        if mol is None:
            raise ValueError("Molecule can't be None.")

        has_hydrogen = False
        for atom in mol.GetAtoms():
            if atom.GetNumExplicitHs() > 0:
                has_hydrogen = True
                break
        return has_hydrogen

    @staticmethod
    def get_inchi_key(mol: Mol) -> str:
        return Chem.MolToInchiKey(mol)

    @staticmethod
    def get_inchi_key_from_smiles(smiles: str) -> str:
        return Chem.MolToInchiKey(
            MoleculeUtils.get_mol_from_smiles(
                MoleculeUtils.get_canonical_smiles(smiles)
            )
        )

    @staticmethod
    def get_molecular_formula(mol: Mol) -> str:
        return CalcMolFormula(mol)

    @staticmethod
    def get_mol_block(mol: Mol, stereo: bool) -> str:
        if stereo:
            hsMol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(hsMol)
            return Chem.MolToMolBlock(hsMol)
        return Chem.MolToMolBlock(mol, True, -1, False)

    @staticmethod
    def get_canonical_smiles_from_mol(mol: Mol, stereo: bool) -> str:
        m = Chem.Mol(mol)
        return Chem.MolToSmiles(standardizer.standardize_mol(m), stereo)

    @staticmethod
    def get_mol_from_smiles(smiles: str) -> Mol:
        return standardizer.standardize_mol(Chem.MolFromSmiles(smiles))

    @staticmethod
    def get_canonical_smiles(smiles: str) -> str:
        return Chem.MolToSmiles(MoleculeUtils.get_mol_from_smiles(smiles), True)

    @staticmethod
    def get_mol_weight(mol: Mol) -> float:
        return Descriptors.MolWt(mol)

    @staticmethod
    def molecules_from_sdf(filepath: str):
        with open(filepath, "rb") as file:
            suppl = Chem.ForwardSDMolSupplier(file)
            for _, rdmol in enumerate(suppl, 1):
                if rdmol:
                    yield rdmol
                else:
                    continue
