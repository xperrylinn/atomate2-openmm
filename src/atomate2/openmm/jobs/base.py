from typing import Optional, Union, Callable
from jobflow import job, Maker, Job
from dataclasses import (
    dataclass,
    field,
)
from pathlib import Path
from pymatgen.io.openmm.sets import OpenMMSet
from pymatgen.io.openmm.inputs import StateInput
from src.atomate2.openmm.schemas.openmm_task_document import (
    OpenMMTaskDocument,
    CalculationInput,
    CalculationOutput,
    PhysicalState,
    TaskDetails,
    StateReports,
    DCDReports,
)
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
from openmm.app.simulation import Simulation
from tempfile import TemporaryDirectory
import copy
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
    # todo: add data keyword argument to specify where to write bigger files like trajectory files
    return job(method, output_schema=OpenMMTaskDocument)


@dataclass
class BaseOpenmmMaker(Maker):
    """
    platform_kwargs : Optional[dict]
        platform and platform_properties passed to OpenMMSet.get_simulation.
        If no arguments provided defaults to CPU.
    dcd_reporter_interval : Optional[int]
        DCD reporter interval. If DCD reporter is not desired, set keyword argument to 0.
    state_reporter_interval : Optional[int]
        State reporter interval. If state reporter is not desired, set keyword argument to 0.
    """
    name: str = "base openmm job"
    platform_kwargs: Optional[dict] = field(default_factory=dict)
    dcd_reporter_interval: Optional[int] = field(default=10000)
    state_reporter_interval: Optional[int] = field(default=1000)
    wrap_dcd: Optional[bool] = False

    @openmm_job
    def make(
            self,
            input_set: OpenMMSet,
            output_dir: Optional[Union[str, Path]] = None,
    ) -> Job:
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
        output_dir : Optional[Union[str, Path]]
            Path to directory for writing state and DCD trajectory files. This could be a temp or
            persistent directory.
        """

        # Define output_dir if as a temporary directory if not provided
        temp_dir = None     # Define potential pointer to temporary directory to keep in scope
        if output_dir is None:
            temp_dir = TemporaryDirectory()
            output_dir = temp_dir.name

        # Setup simulation
        sim = self._setup_base_openmm_task(input_set, output_dir)

        # Run the simulation
        self._run_openmm(sim)

        # Close the simulation
        input_set = self._close_base_openmm_task(input_set, sim.context, output_dir)

        return input_set

    def _setup_base_openmm_task(self, input_set: OpenMMSet, output_dir: Union[str, Path]) -> Simulation:
        """
        Initializes an OpenMM Simulation. Dervied classes define the _run_openmm with simulation task details.

        Parameters
        ----------
        input_set : OpenMMSet
            OpenMM set for initializing an OpenMM Simulation.
        output_dir : Optional[Union[str, Path]]
            Path to directory for writing state and DCD trajectory files. This could be a temp or
            persistent directory.

        Returns
        -------
        sim : Simulation
            OpenMM Simulation from OpenMMSet.

        """

        # Setup compute platform and get a Simulation
        platform_name = self.platform_kwargs.get("platform")
        platform_name = platform_name if platform_name is not None else "CPU"
        platform_props = self.platform_kwargs.get("platform_properties")
        platform = Platform.getPlatformByName(platform_name)

        sim = input_set.get_simulation(
            platform=platform,
            platformProperties=platform_props
        )

        # Add reporters
        if self.dcd_reporter_interval > 0:
            dcd_file_name = os.path.join(output_dir, "trajectory_dcd")
            dcd_reporter = DCDReporter(
                file=dcd_file_name,
                reportInterval=self.dcd_reporter_interval,
            )
            sim.reporters.append(dcd_reporter)
        if self.state_reporter_interval > 0:
            state_file_name = os.path.join(output_dir, "state_xml")
            state_reporter = StateDataReporter(
                file=state_file_name,
                reportInterval=self.state_reporter_interval,
            )
            sim.reporters.append(state_reporter)

        return sim

    def _run_openmm(self, *args, **kwargs) -> TaskDetails:
        """
        Abstract method for holding the logic to be ran by each job.
        """
        raise NotImplementedError("_run_openmm should be implemented by each class that derives from BaseOpenmmMaker.")

    def _close_base_openmm_task(self, input_set: OpenMMSet, context: Context, output_dir: Union[str, Path]):

        # Create an output OpenMMSet for CalculationOutput
        output_set = copy.deepcopy(input_set)   # comment out until - https://github.com/materialsproject/pymatgen/pull/2973/files
        state = StateInput(
            context.getState(
                getPositions=True,
                getVelocities=True,
                enforcePeriodicBox=True,
            )
        )
        output_set[output_set.state_file] = state

        # Grab StateDataReporter and DCDReporter if present on simulation reporters
        state_reports, dcd_reports = None, None
        if self.state_reporter_interval > 0:
            # todo: what happens when state_reporter_interval > 0, but nothing has been
            # reported, for example, Simulation.step was not called. Look at TaskDetails
            # for logic flow?
            state_file_name = os.path.join(output_dir, "state_xml")
            state_reports = StateReports.from_state_file(state_file_name)
        if self.dcd_reporter_interval > 0:
            dcd_file_name = os.path.join(output_dir, "trajectory_dcd")
            dcd_reports = DCDReports(
                location=dcd_file_name,
                report_interval=self.dcd_reporter_interval
            )

        # TODO: create a PDP reporter and report once to create a new PDB file to add to output

        calculation_input = CalculationInput(
            input_set=input_set,
            physical_state=PhysicalState.from_input_set(input_set)
        )
        calculation_output = CalculationOutput(
            output_set=output_set,
            physical_state=PhysicalState.from_input_set(output_set),
            state_reports=state_reports,
            dcd_reports=dcd_reports,
        )

        task_doc = OpenMMTaskDocument(
            input_set=output_set,
            physical_state=PhysicalState.from_input_set(output_set),
            calculation_input=calculation_input,
            calculation_output=calculation_output,
        )

        return task_doc
