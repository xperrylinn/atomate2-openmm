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
    name: str = "OpenMMSolutionGen maker"

    @job
    def make(
            self,
            input_mol_dicts: List[Union[Dict, InputMoleculeSpec]],
            density: Optional[float] = None,
            box: Optional[List[float]] = None,
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
        openmm_sol_gen = OpenMMSolutionGen(**kwargs)
        input_set = openmm_sol_gen.get_input_set(
            input_mol_dicts=input_mol_dicts,
            density=density,
            box=box,
        )

        return input_set


@dataclass
class EnergyMinimizationMaker(BaseOpenmmMaker):
    name: str = "energy minimization"

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Union[str, Path],
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Create a flow that runs in the NPT ensemble.

        Parameters
        ----------
        input_set : OpenMMSet
            A pymatgen structure object.
        platform : Optional[Union[str, openmm.openmm.Platform]]
             platform: the OpenMM platform passed to the Simulation.
        platform_properties : Optional[Dict[str, str]]
            properties of the OpenMM platform that is passed to the simulation.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        platform = Platform.getPlatformByName(platform)
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_properties
        )
        sim.minimizeEnergy()
        new_state_file_path = os.path.join(output_dir, "state_xml")
        sim.saveState(new_state_file_path)
        return input_set


@dataclass
class NPTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Union[str, Path],
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Equilibrate the pressure of a simulation in the NPT ensemble.

        Adds and then removes a openmm.MonteCarlo Barostat to shift the system
        into the NPT ensemble.

        Parameters
        ----------
        input_set : OpenMMSet
            A pymatgen structure object.
        platform : Optional[Union[str, openmm.openmm.Platform]]
             platform: the OpenMM platform passed to the Simulation.
        platform_properties : Optional[Dict[str, str]]
            properties of the OpenMM platform that is passed to the simulation.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """

        platform = Platform.getPlatformByName(platform)
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_properties,
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
        new_state_file_path = os.path.join(output_dir, "state_xml")
        sim.saveState(new_state_file_path)
        return input_set


@dataclass
class NVTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "nvt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Union[str, Path],
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Equilibrate the temperature of a simulation in the NPT ensemble.

        Adds and then removes a openmm.MonteCarlo Barostat to shift the system
        into the NPT ensemble.

        Parameters
        ----------
        input_set : OpenMMSet
            A pymatgen structure object.
        platform : Optional[Union[str, openmm.openmm.Platform]]
             platform: the OpenMM platform passed to the Simulation.
        platform_properties : Optional[Dict[str, str]]
            properties of the OpenMM platform that is passed to the simulation.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        platform = Platform.getPlatformByName(platform)
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_properties,
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
        new_state_file_path = os.path.join(output_dir, "state_xml")
        sim.saveState(new_state_file_path)
        return input_set


@dataclass
class AnnealMaker(BaseOpenmmMaker):

    steps: Tuple[int, int, int] = (500000, 1000000, 500000)
    name: str = "anneal simulation"
    temperatures: Tuple[int, int, int] = (298, 400, 298)
    temp_steps: int = 100

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Union[str, Path],
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
        platform : Optional[Union[str, openmm.openmm.Platform]]
             platform: the OpenMM platform passed to the Simulation.
        platform_properties : Optional[Dict[str, str]]
            properties of the OpenMM platform that is passed to the simulation.

        Returns
        -------
        Job
            A OpenMM job containing one npt run.
        """
        platform = Platform.getPlatformByName(platform)
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_properties,
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
        new_state_file_path = os.path.join(output_dir, "state_xml")
        sim.saveState(new_state_file_path)
        return input_set
