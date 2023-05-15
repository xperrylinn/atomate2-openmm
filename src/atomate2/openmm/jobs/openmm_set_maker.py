from dataclasses import dataclass, field
from tempfile import TemporaryDirectory
from typing import Optional, Union, List, Dict
from pathlib import Path

from jobflow import Maker
from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.schema import InputMoleculeSpec
from pymatgen.io.openmm.sets import OpenMMSet

from src.atomate2.openmm.jobs.base_openmm_maker import openmm_job


@dataclass
class OpenMMSetMaker(Maker):
    """
    Base class for OpenMM set makers.
    """
    name: str = "openmm set maker"
    generator: OpenMMSolutionGen = field(default_factory=OpenMMSolutionGen)

    @openmm_job
    def make(
        self,
        input_mol_dicts: List[Union[Dict, InputMoleculeSpec]],
        density: Optional[float] = None,
        box: Optional[List[float]] = None,
        output_dir: Optional[Union[str, Path]] = None
    ) -> OpenMMSet:
        if output_dir is None:
            temp_dir = TemporaryDirectory()
            output_dir = temp_dir.name
            output_dir = Path(output_dir)
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        output_set = self.generator.get_input_set(input_mol_dicts, density, box)

        output_set.write_input(output_dir)

        return output_set
