from typing import List, Union, Optional, Dict
from dataclasses import dataclass
from jobflow import Maker, job
from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.schema import InputMoleculeSpec


@dataclass
class OpenMMSetFromInputMoleculeSpec(Maker):
    name: str = "OpenMMSetFromInputMoleculeSpec maker"

    @job
    def make(
            self,
            input_mol_dicts: List[Union[Dict, InputMoleculeSpec]],
            density: Optional[float] = None,
            box: Optional[List[float]] = None,
            topology_file: str="topology_pdb",
            state_file: str="state_xml",
            system_file: str="system_xml",
            integrator_file: str="integrator_xml",
            contents_file: str="contents_json",
            **kwargs
    ):
        """
        OpenMMSolutionGen wrapper for generating an OpenMMSet.

        Parameters
        ----------
        input_mol_dicts : List[Union[Dict, InputMoleculeSpec]]
            List of dictionaries or InputMoleculeSpecs
            for passing to OpenMMSolutionGen.
        density : Optional[float]
            Density of simulation, molcules/atoms per cubic centimeter.
            Specify at most one argument between density and box.
        box : Optional[List[float]
            Dimensions of simulation box, in units of todo: ?.
            Specify at most one argument between density and box.
        kwargs : Dict
            Keyword arguments for OpenMMSolutionGen.get_input_set
        Returns
        -------
        Job
            Job for generating an OpenMM input set instance.

        """
        openmm_sol_gen = OpenMMSolutionGen(
            topology_file=topology_file,
            state_file=state_file,
            system_file=system_file,
            integrator_file=integrator_file,
            contents_file=contents_file,
            **kwargs
        )
        input_set = openmm_sol_gen.get_input_set(
            input_mol_dicts=input_mol_dicts,
            density=density,
            box=box,
        )

        return input_set
