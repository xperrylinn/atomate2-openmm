import pathlib
from typing import Union

from pydantic import BaseModel, Field
from pymatgen.io.openmm.sets import OpenMMSet
from pymatgen.io.openmm.schema import SetContents
from atomate2.utils.datetime import datetime_str
# from src.atomate2.openmm import Calculation

from io import StringIO
import numpy as np

from typing import List, Tuple


class StateReports(BaseModel):
    steps: List[int] = Field(None, description="List of steps")
    time: List[float] = Field(None, description="List of times")
    potential_energy: List[float] = Field(None, description="List of potential energies")
    kinetic_energy: List[float] = Field(None, description="List of kinetic energies")
    total_energy: List[float] = Field(None, description="List of total energies")
    temperature: List[float] = Field(None, description="List of temperatures")
    volume: List[float] = Field(None, description="List of volumes")
    density: List[float] = Field(None, description="List of densities")

    @classmethod
    def from_state_file(cls, state_file: Union[str, Path]):
        data = np.loadtxt(state_file, delimiter=',', skiprows=1)

        # Define the order of attributes in the data
        attributes = ['steps', 'potential_energy', 'kinetic_energy', 'total_energy',
                      'temperature', 'volume', 'density']

        # Extract the data columns and set the corresponding class fields
        attributes = {
            attribute: data[:, i].tolist() for i, attribute in enumerate(attributes)
        }

        return StateReports(**attributes)


class DCDReports(BaseModel):
    location: str = Field(None, description="Location of the DCD file")  # this should be a S3 location
    report_interval: int = Field(None, description="Report interval")
    enforce_periodic_box: bool = Field(None, description="Wrap particles or not")

    @classmethod
    def from_dcd_file(cls, dcd_file):
        # TODO: will somehow need to interface with the additional store?
        return


class PhysicalState(BaseModel):
    box_vectors: Tuple[Tuple[float, float, float]] = Field(None, description="Box vectors for the calculation")
    temperature: float = Field(None, description="Temperature for the calculation")
    step_size: float = Field(None, description="Step size for the calculation")
    friction_coefficient: float = Field(None, description="Friction coefficient for the calculation")

    @classmethod
    def from_input_set(cls, input_set):
        integrator = input_set.inputs[input_set.integrator_file].get_integrator()
        state = input_set.inputs[input_set.state_file].get_state()

        vector_array = state.getPeriodicBoxVectors(asNumpy=True)._value
        box_vectors = tuple(tuple(vector) for vector in vector_array)
        temperature = integrator.getTemperature()._value  # kelvin
        step_size = integrator.getStepSize()._value  # picoseconds
        friction_coefficient = integrator.getFriction()._value  # 1/picoseconds

        return PhysicalState(
            box_vectors=box_vectors,
            temperature=temperature,
            step_size=step_size,
            friction_coefficient=friction_coefficient,
        )


class CalculationInput(BaseModel):
    input_set: OpenMMSet = Field(None, description="Input set for the calculation")
    physical_state: PhysicalState = Field(None, description="Physical state for the calculation")
    contents: SetContents = Field(None, description="Contents of the set")

    @classmethod
    def from_input_set(cls, input_set):
        physical_state = PhysicalState.from_input_set(input_set)
        contents = input_set.inputs[input_set.contents_file].contents
        return CalculationInput(
            input_set=input_set,
            physical_state=physical_state,
            contents=contents,
        )


class CalculationOutput(BaseModel):
    input_set: OpenMMSet = Field(None, description="Input set for the calculation")
    physical_state: PhysicalState = Field(None, description="Physical state for the calculation")
    state_reports: StateReports = Field(None, description="State reporter output")
    dcd_reports: DCDReports = Field(None, description="DCD reporter output")

    @classmethod
    def from_directory(cls, output_dir: Union[str, pathlib.Path]):
        # need to write final input_set to output_dir
        # parse state reporter
        # will need to figure out location of dcd_reporter from  additional store, which is global

        # in this approach, we will need to write out the final input_set to the output_dir
        # should we put the OpenMMSet in a sub-directory or just dump it's
        # contents into the output_dir? that will influence this method
        output_dir = pathlib.Path(output_dir)
        input_set = OpenMMSet.from_directory(output_dir)
        return CalculationOutput(
            input_set=input_set,
            physical_state=PhysicalState.from_input_set(input_set),
            # these need to be named consistently when they are written out
            state_reporter=StateReports.from_state_file(output_dir / "???"),
            dcd_reporter=DCDReports.from_dcd_file(output_dir / "???"),
        )


class TaskDetails(BaseModel):
    task_name: str = Field(None, description="Task name")
    task_args: dict = Field(None, description="Task args")
    task_kwargs: dict = Field(None, description="Task kwargs")
    platform: str = Field(None, description="Platform")
    platform_details: dict = Field(None, description="Platform details")
    total_steps: int = Field(None, description="Total steps")


class OpenMMTaskDocument(BaseModel):
    """Definition of the OpenMM task document."""

    output_dir: str = Field(None, description="The directory for this OpenMM task")
    input_dir: OpenMMSet = Field()
    input: CalculationInput = Field(None, description="Input for the calculation")
    output: CalculationOutput = Field(None, description="Output for the calculation")
    details: TaskDetails = Field(None, description="Details about the task")
    elapsed_time: float = Field(None, description="Elapsed time for the calculation")
    # i think we should have this?
    # state: Status = Field(None, description="Status of the calculation")
    last_updated: str = Field(default_factory=datetime_str, description="Timestamp for this task document was last updated")
    completed_at: str = Field(default=None, description="Timestamp for when this task was completed")
    # you are right let's remove this
    # input_set: OpenMMSet = Field(None, description="Final output structure from the task")  # final OpenMM set output obj
    task_label: str = Field(None, description="A description of the task")
    tags: List[str] = Field(None, description="Metadata tags for this task document")
