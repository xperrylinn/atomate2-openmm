from typing import (
    List,
    Union,
    Optional,
    Dict,
    Tuple,
)
from pathlib import Path

from dataclasses import dataclass

from jobflow import Maker, job

from pymatgen.io.openmm.generators import OpenMMSolutionGen
from pymatgen.io.openmm.schema import InputMoleculeSpec
from pymatgen.io.openmm.sets import OpenMMSet
from pymatgen.io.openmm.inputs import StateInput

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
from openmm.app import (
    DCDReporter,
    StateDataReporter,
    PDBReporter,
)
from openmm import Platform

from src.atomate2.openmm.jobs.base import BaseOpenmmMaker

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


@dataclass
class EnergyMinimizationMaker(BaseOpenmmMaker):
    name: str = "energy minimization"

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: str,
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Geometry optimization Job maker.

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
            file=os.path.join(output_dir, "energy_minimation_trajecotry.dcd"),
            reportInterval=10,
        )
        sim.reporters.append(dcd_reporter)

        # Minimize the energy
        sim.minimizeEnergy()

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
            output_dir: str,
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Equilibrate the pressure of a simulation in the NPT ensemble.

        Adds and then removes a openmm.MonteCarloBarostat to shift the system
        into the NPT ensemble.

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
            file=os.path.join(output_dir, "npt_trajecotry.dcd"),
            reportInterval=10,
        )
        sim.reporters.append(dcd_reporter)

        # Add barostat to system
        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        barostat_force_index = system.addForce(
            MonteCarloBarostat(self.pressure * atmosphere, self.temperature * kelvin, 10)
        )

        # Re-init the context afte adding thermostat to System
        context.reinitialize(preserveState=True)

        # Run the simulation
        sim.step(self.steps)

        # Remove thermostat and update context
        system.removeForce(barostat_force_index)
        context.reinitialize(preserveState=True)

        # output_set = copy.deepcopy(input_set)
        state = StateInput(
            context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=True,
            )
        )

        input_set[input_set.state_file] = state
        return input_set


@dataclass
class NVTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "nvt simulation"
    temperature: float = 298
    frequency: int = 10

    def _run_openmm(self, sim):    # parameters should be NVTMaker attrs
        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        thermostat_force_index = system.addForce(
            AndersenThermostat(self.temperature * kelvin, pico * second * self.frequency)
        )

        # Re-init the context afte adding thermostat to System
        context.reinitialize(preserveState=True)

        # Run the simulation
        sim.step(self.steps)

        # Remove thermostat and update context
        system.removeForce(thermostat_force_index)
        context.reinitialize(preserveState=True)

    @job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: str,
            platform: Optional[Union[str, Platform]] = "CPU",
            platform_properties: Optional[Dict[str, str]] = None,
    ):
        """
        Equilibrate the temperature of a simulation in the NVT ensemble.

        Adds and then removes a openmm.AndersenThermostat Thermostat to shift the system
        into the NVT ensemble.

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
            A OpenMM job containing one nvt run.
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
            file=os.path.join(output_dir, "nvt_trajecotry.dcd"),
            reportInterval=10,
        )
        sim.reporters.append(dcd_reporter)


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


@dataclass
class AnnealMaker(BaseOpenmmMaker):
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
