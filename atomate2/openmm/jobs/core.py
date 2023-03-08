
from typing import List, Union, Optional, Dict
from pathlib import Path

from dataclasses import dataclass

from jobflow import Maker

from pymatgen.io.openmm.schema import InputMoleculeSpec
from pymatgen.io.openmm.sets import OpenMMSet

from atomate2.openmm.jobs.base import BaseOpenmmMaker


@dataclass
class InputMaker(Maker):
    name: str = "input maker"

    def make(
        self,
        input_mol_dicts: List[Union[Dict, InputMoleculeSpec]],
        density: Optional[float] = None,
        box: Optional[List[float]] = None,
        prev_dir: str | Path | None = None
    ):
        """
        Create a job that generates an OpenMM input set.

        Parameters
        ----------
        input_set : .OpenMMSet
            A pymatgen structure object.
        density : float
        box : list
        prev_dir : str or Path or None
            A previous OpenMM simulation directory to link files from.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        return


@dataclass
class EnergyMinimizationMaker(BaseOpenmmMaker):
    name: str = "energy minimization"

    def make(self, input_set: OpenMMSet, prev_dir: str | Path | None = None):
        """
        Create a flow that runs in the NPT ensemble.

        Parameters
        ----------
        input_set : .OpenMMSet
            A pymatgen structure object.
        prev_dir : str or Path or None
            A previous OpenMM simulation directory to link files from.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        return



@dataclass
class NPTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    def make(self, input_set: OpenMMSet, prev_dir: str | Path | None = None):
        """
        Create a flow that runs in the NPT ensemble.

        Parameters
        ----------
        input_set : .OpenMMSet
            A pymatgen structure object.
        prev_dir : str or Path or None
            A previous OpenMM simulation directory to link files from.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        return


@dataclass
class NVTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    def make(self, input_set: OpenMMSet, prev_dir: str | Path | None = None):
        """
        Create a flow that runs in the NPT ensemble.

        Parameters
        ----------
        input_set : .OpenMMSet
            A pymatgen structure object.
        prev_dir : str or Path or None
            A previous OpenMM simulation directory to link files from.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        return


@dataclass
class AnnealMaker(BaseOpenmmMaker):

    steps: tuple[int, int, int] = (500000, 1000000, 500000)
    name: str = "npt simulation"
    temperature: tuple[int, int, int] = (298, 400, 298)

    def make(self, input_set: OpenMMSet, prev_dir: str | Path | None = None):
        """
        Create a flow that runs in the NPT ensemble.

        Parameters
        ----------
        input_set : .OpenMMSet
            A pymatgen structure object.
        prev_dir : str or Path or None
            A previous OpenMM simulation directory to link files from.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        return
