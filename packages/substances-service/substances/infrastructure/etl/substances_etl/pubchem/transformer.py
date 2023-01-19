import concurrent.futures
import gzip
import logging
import pickle  # nosec
from dataclasses import dataclass
from itertools import islice
from typing import Any, List, Union

import jsonpickle  # type: ignore
from rdkit.Chem import Mol  # type: ignore
from sciences.chemistry.molecules import MoleculeUtils
from shared.utils.path import PathUtils
from substances_core.domain.susbtances.objects import (
    CanonicalSmiles,
    Formula,
    InchiKey,
    IupacName,
    SubstanceName,
)
from substances_core.domain.susbtances.substance import Substance
from substances_core.domain.syncs.entities.item import SyncItem
from substances_core.domain.syncs.services.transformer import SyncTransformer

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class PubchemSubstanceTransformerOptions:
    """PubchemTransformerOptions"""

    data_dir: str
    """Path to the directory used to store extracted data."""

    trf_extension: str = "pickle"
    """"""

    id_name: str = "PUBCHEM_COMPOUND_CID"
    """Key used to identify the compound ID insides SDF."""

    smiles_name: str = "PUBCHEM_OPENEYE_ISO_SMILES"

    mol_name: str = "Molecule"

    chunk_size: int = 100


class PuchemSubstanceTransformer(SyncTransformer[Substance]):
    """PuchemSubstanceTransformer"""

    def __init__(self, options: PubchemSubstanceTransformerOptions) -> None:
        self._options = options
        # self._warnings = []

    def transform(self, item: SyncItem, inmemory: bool = False) -> Union[str, bytes]:
        try:
            entities = self.__parse_parallel(
                PathUtils.join_path(self._options.data_dir, item.token),
                self._options.chunk_size,
            )

            compressed = gzip.compress(pickle.dumps(jsonpickle.encode(entities)))
            if inmemory:
                return compressed
            else:
                pickle_file = f"{item.token}.{self._options.trf_extension}"
                output_path = PathUtils.join_path(self._options.data_dir, pickle_file)
                with open(output_path, "wb") as f_out:
                    f_out.write(compressed)
                return pickle_file
        except Exception as err:
            LOGGER.error(err)
            raise

    def parse(self, data: Any, chunk_start=0, chunk_size=100) -> List[Substance]:
        print(f"Parse chunk start: {chunk_start}-{chunk_start + chunk_size}")
        substances: List[Substance] = []
        generator = MoleculeUtils.molecules_from_sdf(data)

        for molecule in islice(generator, chunk_start, chunk_start + chunk_size):
            try:
                if MoleculeUtils.has_explicit_hygrogen(molecule):
                    substance = Substance(
                        inchi_key=InchiKey(MoleculeUtils.get_inchi_key(molecule)),
                        canonical_smiles=CanonicalSmiles(
                            self.__extract_canonical_smiles(molecule, True)
                        ),
                        name=SubstanceName(self.__extract_prefered_name(molecule)),
                        iupac_name=IupacName(self.__extract_iupac_name(molecule)),
                        formula=Formula(MoleculeUtils.get_molecular_formula(molecule)),
                        mol_weight=self.__extract_mol_weight(molecule),
                        mass=self.__extract_mass(molecule),
                        mol3d=MoleculeUtils.get_mol_block(molecule, True),
                    )
                else:
                    substance = Substance(
                        inchi_key=InchiKey(MoleculeUtils.get_inchi_key(molecule)),
                        canonical_smiles=CanonicalSmiles(
                            self.__extract_canonical_smiles(molecule, False)
                        ),
                        name=SubstanceName(self.__extract_prefered_name(molecule)),
                        iupac_name=IupacName(self.__extract_iupac_name(molecule)),
                        formula=Formula(MoleculeUtils.get_molecular_formula(molecule)),
                        mol_weight=self.__extract_mol_weight(molecule),
                        mass=self.__extract_mass(molecule),
                        mol2d=MoleculeUtils.get_mol_block(molecule, False),
                    )

                substances.append(substance)
            except Exception as err:
                LOGGER.error(err)
                print(f"Molecule transformation error: {dict(molecule)}")

        print(f"Parse chunk done: {chunk_start}-{chunk_start + chunk_size}")

        return substances

    def __parse_parallel(self, data, chunk_size=100):
        generator = MoleculeUtils.molecules_from_sdf(data)
        item_count: int = sum(1 for _ in generator)
        print(f"Items to parse: {item_count}")
        try:
            substances = []
            # TODO: import check for choose correct concurrent method - > pool is not available on AWS Lambda
            # with multiprocessing.Pool(multiprocessing.cpu_count()) as pool:
            #     jobs = []
            #     for i in range(0, item_count, chunk_size):
            #         jobs.append(pool.apply_async(self.parse, (data, i, chunk_size)))
            #     for job in jobs:
            #         result = job.get()
            #         substances = list(itertools.chain(substances, result))
            with concurrent.futures.ThreadPoolExecutor(max_workers=32) as executor:
                futures = [
                    executor.submit(self.parse, data, i, chunk_size)
                    for i in range(0, item_count, chunk_size)
                ]
                for future in concurrent.futures.as_completed(futures):
                    substances = [*substances, *future.result()]
        except Exception as err:
            LOGGER.error(err)
            raise
        finally:
            return substances

    def __extract_canonical_smiles(self, molecule: Mol, stereo: bool) -> str:
        return MoleculeUtils.get_canonical_smiles_from_mol(molecule, stereo)

    def __extract_mol_weight(self, molecule: Mol) -> float:
        mol_weight: float = molecule.GetProp("PUBCHEM_MOLECULAR_WEIGHT")
        return mol_weight

    def __extract_prefered_name(self, molecule: Mol) -> str:
        name: str = molecule.GetProp("PUBCHEM_IUPAC_NAME")
        return name

    def __extract_iupac_name(self, molecule: Mol) -> str:
        name: str = molecule.GetProp("PUBCHEM_IUPAC_NAME")
        return name

    def __extract_mass(self, molecule: Mol) -> float:
        mass: float = molecule.GetProp("PUBCHEM_EXACT_MASS")
        return float(mass)
