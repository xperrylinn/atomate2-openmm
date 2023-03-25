from typing import List, Union, Optional, Dict, Tuple
from pathlib import Path

from dataclasses import dataclass

from jobflow import Maker, job

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.schema import InputMoleculeSpec
from pymatgen.io.openmm.sets import OpenMMSet

from openmm.openmm import (
    MonteCarloBarostat,
    AndersenThermostat,
)
from openmm.unit import (
    kelvin,
    atmosphere,
    pico,
    second
)
from openmm import Platform

from src.atomate2.openmm.jobs.base import BaseOpenmmMaker

import numpy as np


@dataclass
class InputMaker(Maker):
    name: str = "input maker"

    @job
    def make(
        self,
        input_mol_dicts: List[Union[Dict, InputMoleculeSpec]],
        density: Optional[float] = None,
        box: Optional[List[float]] = None,
        prev_dir: Optional[Union[str, Path]] = None
    ):
        """
        Create a job for generating an OpenMMSet instance

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
        prev_dir : Optional[Union[str, Path]]
            Previous OpenMM simulation directory to link files from.
        Returns
        -------
        Job
            Job for generating an OpenMM input set instance.

        """
        input_set = OpenMMSolutionGen().get_input_set(
            input_mol_dicts=input_mol_dicts,
            density=density,
            box=box,
        )
        return input_set


@dataclass
class EnergyMinimizationMaker(BaseOpenmmMaker):
    name: str = "energy minimization"

    @job
    def make(self, input_set: OpenMMSet, prev_dir: Optional[Union[str, Path]] = None):
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
        platform = Platform.getPlatformByName("OpenCL")
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties={"DeviceIndex": str(0)},
        )
        sim.minimizeEnergy()
        return input_set


@dataclass
class NPTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    @job
    def make(self, input_set: OpenMMSet, prev_dir: Optional[Union[str, Path]] = None):
        """
        Equilibrate the pressure of a simulation in the NPT ensemble.

        Adds and then removes a openmm.MonteCarlo Barostat to shift the system
        into the NPT ensemble.

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

        platform = Platform.getPlatformByName("OpenCL")
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties={"DeviceIndex": str(0)},
        )

        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        barostat_force_index = system.addForce(
            MonteCarloBarostat(self.pressure * atmosphere, self.temperature * kelvin, 10)
        )
        context.reinitialize(preserveState=True)
        sim.step(self.steps)
        system.removeForce(barostat_force_index)
        context.reinitialize(preserveState=True)

        return input_set


@dataclass
class NVTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "nvt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    @job
    def make(self, input_set: OpenMMSet, prev_dir: Optional[Union[str, Path]] = None):
        """
        Equilibrate the temperature of a simulation in the NPT ensemble.

        Adds and then removes a openmm.MonteCarlo Barostat to shift the system
        into the NPT ensemble.

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
        platform = Platform.getPlatformByName("OpenCL")
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties={"DeviceIndex": str(0)},
        )

        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        thermostat_force_index = system.addForce(
            AndersenThermostat(self.temperature * kelvin, self.frequency * pico * second)
        )
        context.reinitialize(preserveState=True)
        sim.step(self.steps)
        system.removeForce(thermostat_force_index)
        context.reinitialize(preserveState=True)

        return input_set


@dataclass
class AnnealMaker(BaseOpenmmMaker):

    steps: Tuple[int, int, int] = (500000, 1000000, 500000)
    name: str = "anneal simulation"
    temperatures: Tuple[int, int, int] = (298, 400, 298)
    temp_steps: int = 100

    @job
    def make(self, input_set: OpenMMSet, prev_dir: Optional[Union[str, Path]] = None):
        """
        Anneal the simulation at the specified temperature.

        Annealing takes place in 3 stages, heating, holding, and cooling. The three
        elements of steps specify the length of each stage. After heating, and holding,
        the system will cool to its original temperature.

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
        platform = Platform.getPlatformByName("OpenCL")
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties={"DeviceIndex": str(0)},
        )

        # TODO: timing is currently bugged and not propagating long enough, should be fixed
        assert len(self.steps) == 3, ""
        integrator = sim.context.getIntegrator()
        old_temperature = integrator.getTemperature()
        temp_step_size = abs(self.temperature * kelvin - old_temperature) / self.temp_steps

        for temp in self.temperatures:
            integrator.setTemperature(temp * kelvin)
            sim.step(self.steps[0] // self.temp_steps)

        sim.step(self.steps[1])

        for temp in self.temperatures:
            integrator.setTemperature(temp * kelvin)
            sim.step(self.steps[2] // self.temp_steps)

        return input_set
