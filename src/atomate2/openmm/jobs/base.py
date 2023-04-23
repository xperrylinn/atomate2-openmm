from typing import Optional, Union, Callable
from jobflow import job, Maker
from dataclasses import dataclass
from pathlib import Path
from pymatgen.io.openmm.sets import OpenMMSet
from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
from openmm import (
    Platform,
    Context,
)
from typing import (
    Union,
    Optional,
    Dict,
)
from openmm.app import (
    DCDReporter,
    StateDataReporter,
    PDBReporter,
)
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
import os


def openmm_job(method: Callable):
    """

    Parameters
    ----------
    method : callable
        A BaseOpenmmMaker.make method. This should not be specified directly and is
        implied by the decorator.

    Returns
    -------
    callable
        A decorated version of the make function that will generate OpenMM jobs.
    """
    return job(method, output_schema=OpenMMTaskDocument)


@dataclass
class BaseOpenmmMaker(Maker):
    """
    todo: doc string
    """
    name: str = "base openmm job"
    get_simulation_kwargs: Optional[dict] = None
    dcd_reporter_kwargs: Optional[dict] = None
    pdb_reporter_kwargs: Optional[dict] = None
    state_reporter_kwargs: Optional[dict] = None
    pressure_barostat_kwargs: Optional[dict] = None
    temperature_barostat_kwargs: Optional[dict] = None

    @openmm_job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Optional[Union[str, Path]] = None,
    ):
        """
        Run an Openmm calculation.

        this should mirror the structure of BaseVaspMaker.make
        it should include everything needed to setup and teardown and openmm
        simulation but it should offload the details to its children
        (e.g. NPTMaker, NVTMaker, etc.)

        Parameters
        ----------
        input_set : OpenMMSet
            pymatgen.io.openmm OpenMMSet object instance.
        prev_dir : str or Path or None
            Previous OpenMM calculation directory to copy output files from.
        """

        # Setup simulation
        sim, foces_indices = self._setup_base_openmm_task(input_set, output_dir)

        # Run the simulation
        sim.step(self.steps)

        # Close the simulation
        input_set = self._close_base_openmm_task()

        return input_set

    def _setup_base_openmm_task(self, input_set: OpenMMSet, output_dir: Optional[Union[str, Path]] = None):

        # Setup compute platform and get a Simulation
        platform_name = self.get_simulation_kwargs.get("platform")
        platform = Platform.getPlatformByName(platform_name if platform_name else "CPU")
        platform_props = self.get_simulation_kwargs.get("platform_properties")
        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_props
        )

        # Add reporters
        if self.dcd_reporter_kwargs:
            self.dcd_reporter_kwargs["file"] = os.path.join(output_dir, self.dcd_reporter_kwargs["file"])
            dcd_reporter = DCDReporter(**self.dcd_reporter_kwargs)
            sim.reporters.append(dcd_reporter)
        if self.state_reporter_kwargs:
            self.state_reporter_kwargs["file"] = os.path.join(output_dir, self.state_reporter_kwargs["file"])
            state_reporter = StateDataReporter(**self.state_reporter_kwargs)
            sim.reporters.append(state_reporter)
        if self.pdb_reporter_kwargs:
            self.pdb_reporter_kwargs["file"] = os.path.join(output_dir, self.pdb_reporter_kwargs["file"])
            pdb_reporter = DCDReporter(**self.pdb_reporter_kwargs)
            sim.reporters.append(pdb_reporter)

        # Add barotats
        force_indices = list()
        context = sim.context
        system = context.getSystem()
        if self.temperature_barostat_kwargs:
            thermostat_force_index = system.addForce(
                AndersenThermostat(self.temperature * kelvin, pico * second * self.frequency)
            )
            context.reinitialize(preserveState=True)
            force_indices.append(thermostat_force_index)
        if self.pressure_barostat_kwargs:
            assert (
                system.usesPeriodicBoundaryConditions()
            ), "system must use periodic boundary conditions for pressure equilibration."
            barostat_force_index = system.addForce(
                MonteCarloBarostat(self.pressure * atmosphere, self.temperature * kelvin, 10)
            )
            context.reinitialize(preserveState=True)
            force_indices.append(barostat_force_index)

        return sim, force_indices

    def _close_base_openmm_task(self, context: Context, force_indices: Optional[list]):
        # output_set = copy.deepcopy(input_set)
        state = StateInput(
            context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=True,
            )
        )

        input_set[input_set.state_file] = state
        task_doc =
        return input_set
