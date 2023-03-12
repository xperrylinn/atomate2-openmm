
from pydantic import BaseModel, Field
from pymatgen.io.openmm.sets import OpenMMSet
from atomate2.utils.datetime import datetime_str
from atomate2.openmm.calculation import Calculation

from typing import List
from datetime import datetime


class OpenMMTaskDocument(BaseModel):
    """Definition of the OpenMM task document."""

    dir_name: str = Field(None, description="The directory for this VASP task")
    last_updated: str = Field(
        default_factory=datetime_str,
        description="Timestamp for this task document was last updated",
    )
    completed_at: str = Field(
        None, description="Timestamp for when this task was completed"
    )
    structure: OpenMMSet = Field(
        None, description="Final output structure from the task"
    )

    task_label: str = Field(None, description="A description of the task")
    tags: List[str] = Field(None, description="Metadata tags for this task document")
    calcs_reversed: List[Calculation] = Field(
        None, description="The inputs and outputs for all VASP runs in this task."
    )
