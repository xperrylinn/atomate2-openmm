from jobflow import Maker
from dataclasses import dataclass
from pathlib import Path
from typing import Union
from pymatgen.io.openmm.sets import OpenMMSet
from jobflow import job


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
