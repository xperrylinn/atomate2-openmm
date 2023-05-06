from pydantic import BaseModel, Field
from atomate2.utils.datetime import datetime_str
from src.atomate2.openmm.schemas.calculation_input import CalculationInput
from src.atomate2.openmm.schemas.calculation_output import CalculationOutput
from src.atomate2.openmm.schemas.task_details import TaskDetails
from typing import List, Tuple


class OpenMMTaskDocument(BaseModel):
    """Definition of the OpenMM task document."""

    output_dir: str = Field(None, description="The directory for this OpenMM task")
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
