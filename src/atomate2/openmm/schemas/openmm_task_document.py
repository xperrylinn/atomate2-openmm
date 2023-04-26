from pydantic import BaseModel, Field
from pymatgen.io.openmm.sets import OpenMMSet
from atomate2.utils.datetime import datetime_str
# from src.atomate2.openmm import Calculation

from typing import List


class OpenMMTaskDocument(BaseModel):
    """Definition of the OpenMM task document."""

    output_dir: str = Field(None, description="The directory for this OpenMM task")
    input_dir: OpenMMSet = Field()
    # todo: input summary
    # todo: output summary
    last_updated: str = Field(default_factory=datetime_str, description="Timestamp for this task document was last updated")
    completed_at: str = Field(default=None, description="Timestamp for when this task was completed")
    input_set: OpenMMSet = Field(None, description="Final output structure from the task") # final OpenMM set output obj
    task_label: str = Field(None, description="A description of the task")
    tags: List[str] = Field(None, description="Metadata tags for this task document")
