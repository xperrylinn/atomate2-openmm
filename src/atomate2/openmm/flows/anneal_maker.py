from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from typing import Union, Optional, Dict, Tuple
from pymatgen.io.openmm.sets import OpenMMSet
from openmm import Platform
from dataclasses import dataclass
from jobflow import job
from openmm.app import DCDReporter
from openmm.unit import kelvin
import os
import numpy as np
from pymatgen.io.openmm.inputs import StateInput

"""
    TODO: Refactor into a Flow
"""

@dataclass
class AnnealMaker(BaseOpenMMMaker):
    """
    steps : Union[Tuple[int, int, int], int]

    """

    steps: Union[Tuple[int, int, int], int] = (500000, 1000000, 500000)
    name: str = "anneal simulation"
    temperatures: Union[Tuple[int, int, int], int] = (298, 400, 298)
    temp_steps: int = 100

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
        # Set compute platform
        platform = Platform.getPlatformByName(platform)

        # Get a fresh Simulation from OpenMMSet input argument
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_properties
        )

        # Add DCD reporter
        dcd_reporter = DCDReporter(
            file=os.path.join(output_dir, "anneal_trajecotry.dcd"),
            reportInterval=10,
        )
        sim.reporters.append(dcd_reporter)

        context = sim.context
        integrator = sim.context.getIntegrator()
        old_temperature = integrator.getTemperature()

        # Type checking to assemble annealing step count and temperatures for heating, holding, cooling
        if isinstance(self.temperatures, int):
            anneal_temps = [old_temperature, self.temperatures, old_temperature]
        else:
            anneal_temps = self.temperatures
        if isinstance(self.temp_steps, int):
            anneal_steps = [self.temp_steps] * 3
        else:
            anneal_steps = self.temp_steps

        # Heating temperature
        temp_step_size = abs(anneal_temps[0] * kelvin - old_temperature) / anneal_steps[0]
        for temp in np.arange(
                old_temperature + temp_step_size,
                anneal_temps[0] * kelvin + temp_step_size,
                temp_step_size,
        ):
            integrator.setTemperature(temp * kelvin)
            sim.step(self.steps[0] // anneal_steps[0])

        # Holding
        temp_step_size = anneal_steps[1]
        sim.step(self.steps[1])

        # Cooling temperature
        temp_step_size = abs(anneal_temps[1] * kelvin - anneal_temps[2] * kelvin) / anneal_steps[2]
        for temp in np.arange(
                anneal_temps[1] * kelvin - temp_step_size,
                anneal_temps[2] * kelvin - temp_step_size + temp_step_size,
                -1 * temp_step_size,
        ):
            integrator.setTemperature(temp * kelvin)
            sim.step(self.steps[2] // anneal_steps[2])

        # Get a fresh state object and update the OpenMMSet
        state = StateInput(
            sim.context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=True,
            )
        )
        input_set[input_set.state_file] = state

        return input_set
