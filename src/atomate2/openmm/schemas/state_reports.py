from pydantic import BaseModel, Field
from typing import List, Union
import numpy as np
import pathlib
import logging


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
    def from_state_file(cls, state_file: Union[str, pathlib.Path]):
        data = np.loadtxt(state_file, delimiter=',', skiprows=1)

        if len(data) == 0:
            logging.warning(f"The loaded state file: {state_file}, was empty")
        else:
            # Define the order of attributes in the data
            attributes = ['steps', 'potential_energy', 'kinetic_energy', 'total_energy',
                          'temperature', 'volume', 'density']

            # Extract the data columns and set the corresponding class fields
            attributes = {
                attribute: data[:, i].tolist() for i, attribute in enumerate(attributes)
            }

        return StateReports(**attributes)
