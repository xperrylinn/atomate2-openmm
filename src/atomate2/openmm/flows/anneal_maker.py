from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from typing import Union, Optional, Dict, Tuple
from pymatgen.io.openmm.sets import OpenMMSet
from openmm import Platform
from dataclasses import dataclass
from jobflow import job, Flow
from openmm.app import DCDReporter
from openmm.unit import kelvin
import os
import numpy as np
from pydantic import Field
from pymatgen.io.openmm.inputs import StateInput

from src.atomate2.openmm.jobs.npt_maker import NPTMaker
from src.atomate2.openmm.jobs.temperature_maker import TempChangeMaker

"""
    TODO: Refactor into a Flow
"""

@dataclass
class AnnealMaker(BaseOpenMMMaker):
    """
    steps : Union[Tuple[int, int, int], int]

    """

    name: str = "anneal"
    raise_temp_maker: TempChangeMaker = Field(default_factory=lambda: TempChangeMaker(final_temp=400))
    npt_maker: NPTMaker = Field(default_factory=lambda: NPTMaker())
    lower_temp_maker: TempChangeMaker = Field(default_factory=lambda: TempChangeMaker())

    @staticmethod
    def from_temps_and_steps(
            anneal_temp: int = 400,
            final_temp: int = 298,
            steps: Union[int, Tuple[int, int, int]] = 1500000,
            temp_steps: Union[int, Tuple[int, int, int]] = 100,
            names: Tuple[str, str, str] = ("raise temp", "hold temp", "lower temp")
    ):
        if isinstance(steps, int):
            steps = (steps // 3, steps // 3, steps - 2 * (steps // 3))
        if isinstance(temp_steps, int):
            temp_steps = (temp_steps, temp_steps, temp_steps)

        raise_temp_maker = TempChangeMaker(
            steps=steps[0],
            name=names[0],
            final_temp=anneal_temp,
            temp_steps=temp_steps[0]
        )
        npt_maker = NPTMaker(
            steps=steps[1],
            name=names[1],
        )
        lower_temp_maker = TempChangeMaker(
            steps=steps[2],
            name=names[2],
            final_temp=final_temp,
            temp_steps=temp_steps[2]
        )
        return AnnealMaker(
            raise_temp_maker=raise_temp_maker,
            npt_maker=npt_maker,
            lower_temp_maker=lower_temp_maker,
        )


    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: str,
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Anneal the simulation at the specified temperature.

        Annealing takes place in 3 stages, heating, holding, and cooling. The three
        elements of steps specify the length of each stage. After heating, and holding,
        the system will cool to its original temperature.

        Parameters
        ----------
        input_set : OpenMMSet
            A pymatgen structure object.
        output_dir : str
            Directory to write reporter files to.
        platform : Optional[Union[str, openmm.openmm.Platform]]
             platform: the OpenMM platform passed to the Simulation.
        platform_properties : Optional[Dict[str, str]]
            properties of the OpenMM platform that is passed to the simulation.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        raise_temp_job = self.raise_temp_maker.make(
            input_set=input_set,
            output_dir=output_dir,
        )
        npt_job = self.npt_maker.make(
            input_set=input_set,
            output_dir=output_dir,
        )
        lower_temp_job = self.lower_temp_maker.make(
            input_set=input_set,
            output_dir=output_dir,
        )

        flow = Flow([raise_temp_job, npt_job, lower_temp_job])

        return flow
