from chembl_structure_pipeline import standardizer
from rdkit import Chem

# SAMPLE:
# METACETAMOL
# CANONICAL (INTERNAL) =
# https://pubchem.ncbi.nlm.nih.gov/compound/12124 - CC(=O)NC1=CC(=CC=C1)O
# https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL9419/ - CC(=O)Nc1cccc(O)c1
# https://drugs.ncats.io/drug/V942ZCN81H - CC(=Nc1cccc(c1)O)O (Chemeo not found)
# https://medkoo.com/products/26852 - CC(NC1=CC=CC(O)=C1)=O


def test_standardizer():
    smis = ["CC(=O)NC1=CC(=CC=C1)O", "CC(=O)Nc1cccc(O)c1"]

    cans = [
        Chem.MolToSmiles(standardizer.standardize_mol(Chem.MolFromSmiles(smi)))
        for smi in smis
    ]
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolToSmiles
    # https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.MolFromSmiles
    blocks = [Chem.MolToMolBlock(Chem.MolFromSmiles(can)) for can in cans]

    td_molblock1 = standardizer.standardize_molblock(blocks[0])
    td_molblock2 = standardizer.standardize_molblock(blocks[1])

    assert td_molblock1 == td_molblock2


def test_canonical_smiles():
    smis = [
        "CC(=O)NC1=CC(=CC=C1)O",
        "CC(=O)Nc1cccc(O)c1",
        "CC(=Nc1cccc(c1)O)O",
        "CC(NC1=CC=CC(O)=C1)=O",
    ]

    cans = [
        Chem.MolToSmiles(standardizer.standardize_mol(Chem.MolFromSmiles(smi)))
        for smi in smis
    ]
    assert cans[0] == cans[1]
