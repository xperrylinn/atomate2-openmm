from typing import List, Union, Optional, Dict, Tuple
from pathlib import Path
from dataclasses import dataclass
from jobflow import Maker, job
from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.schema import InputMoleculeSpec
from pymatgen.io.openmm.sets import OpenMMSet
from pymatgen.io.openmm.inputs import StateInput
from openmm.unit import kelvin, atmosphere, pico, second
from openmm.app import DCDReporter
from openmm import Platform
from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenmmMaker
import numpy as np
import os


@dataclass
class OpenMMSetFromDirectory(Maker):
    name: str = "OpenMMSetFromDirectory maker"

    @job
    def make(
            self,
            input_dir: Union[str, Path],
            **kwargs
    ):
        """
        OpenMMSet.from_directory wrapper for generating an OpenMMSet from a directory.

        Parameters
        ----------
        input_dir : Union[str, Path]
            Path to directory containing topology, system, integrator, and state files.
            See OpenMMSet.from_directory for default file names.
        kwargs : Dict
            Keyword arguments for OpenMMSet.from_directory
        Returns
        -------
        Job
            Job for generating an OpenMM input set instance.

        """
        input_set = OpenMMSet.from_directory(
            directory=input_dir,
            topology_file="topology_pdb",
            state_file="state_xml",
            system_file="system_xml",
            integrator_file="integrator_xml",
        )
        return input_set


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
