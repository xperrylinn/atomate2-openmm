from pydantic import BaseModel, Field
from pymatgen.io.openmm.sets import OpenMMSet
from atomate2.utils.datetime import datetime_str
# from src.atomate2.openmm import Calculation

from typing import List


class StateReporterOutput(BaseModel):
    steps: List[int] = Field(None, description="List of steps")
    time: List[float] = Field(None, description="List of times")
    potential_energy: List[float] = Field(None, description="List of potential energies")
    kinetic_energy: List[float] = Field(None, description="List of kinetic energies")
    total_energy: List[float] = Field(None, description="List of total energies")
    temperature: List[float] = Field(None, description="List of temperatures")
    volume: List[float] = Field(None, description="List of volumes")
    density: List[float] = Field(None, description="List of densities")

    @classmethod
    def from_state_file(cls, state_str):
        # TODO: write this once we have a state.txt output to test on
        pass
    pass


class DCDReporterOutput(BaseModel):
    location: str = Field(None, description="Location of the DCD file")  # this should be a S3 location
    report_interval: int = Field(None, description="Report interval")
    enforce_periodic_box: bool = Field(None, description="Wrap particles or not")
    pass


class OpenMMSetContents(BaseModel):
    # this should be in pymatgen.io.openmm
    molecule_specs: list = Field(None, description="Components for the calculation")
    force_fields: List[str] = Field(None, description="Force field for the calculation")
    partial_charge_methods: List[str] = Field(None, description="Partial charge method for the calculation")
    atom_types: List[int] = Field(None, description="Atom types for the calculation")
    atom_resnames: List[str] = Field(None, description="Atom resnames for the calculation")

    @classmethod
    def from_input_set(cls, input_set):
        # get settings from settings
        pass


class PhysicalState(BaseModel):
    box_vectors: list = Field(None, description="Box vectors for the calculation")
    temperature: float = Field(None, description="Temperature for the calculation")
    step_size: float = Field(None, description="Step size for the calculation")
    friction_coefficient: float = Field(None, description="Friction coefficient for the calculation")

    @classmethod
    def from_input_set(cls, input_set):
        # get temp and stuff from integrator
        # get box from state
        pass


class CalculationInput(BaseModel):
    input_set: OpenMMSet = Field(None, description="Input set for the calculation")
    physical_state: PhysicalState = Field(None, description="Physical state for the calculation")
    contents: OpenMMSetContents = Field(None, description="Contents of the set")

    @classmethod
    def from_input_set(cls, input_set):
        # get settings from settings
        # get temp and stuff from integrator
        # get box from state
        pass
    pass


class CalculationOutput(BaseModel):
    input_set: OpenMMSet = Field(None, description="Input set for the calculation")
    physical_state: PhysicalState = Field(None, description="Physical state for the calculation")
    state_reporter: StateReporterOutput = Field(None, description="State reporter output")
    dcd_reporter: DCDReporterOutput = Field(None, description="DCD reporter output")

    @classmethod
    def from_output_dir(cls, output_dir):
        # need to write final input_set to output_dir
        # parse state reporter
        # will need to figure out location of dcd_reporter from
        # the additional stores, which is a global setting
        pass
    pass


class TaskDetails(BaseModel):
    task_name = Field(None, description="Task name")
    task_args = Field(None, description="Task args")
    task_kwargs = Field(None, description="Task kwargs")
    platform = Field(None, description="Platform")
    platform_details = Field(None, description="Platform details")
    total_steps = Field(None, description="Total steps")


class OpenMMTaskDocument(BaseModel):
    """Definition of the OpenMM task document."""

    output_dir: str = Field(None, description="The directory for this OpenMM task")
    input_dir: OpenMMSet = Field()
    input: CalculationInput = Field(None, description="Input for the calculation")
    output: CalculationOutput = Field(None, description="Output for the calculation")
    details: TaskDetails = Field(None, description="Details about the task")
    elapsed_time: float = Field(None, description="Elapsed time for the calculation")
    # state: Status = Field(None, description="Status of the calculation")
    last_updated: str = Field(default_factory=datetime_str, description="Timestamp for this task document was last updated")
    completed_at: str = Field(default=None, description="Timestamp for when this task was completed")
    # you are right let's remove this
    # input_set: OpenMMSet = Field(None, description="Final output structure from the task")  # final OpenMM set output obj
    task_label: str = Field(None, description="A description of the task")
    tags: List[str] = Field(None, description="Metadata tags for this task document")
