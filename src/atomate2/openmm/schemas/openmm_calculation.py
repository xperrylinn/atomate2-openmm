from pydantic import BaseModel


class OpenMMCalculation(BaseModel):
    """
    Full OpenMM calculation inputs and outputs.

    This should contain several core piece of information.

    1. The initial InputSet and an InputSet representing the state at the end of
    the calculation.

    2. A link to the trajectory file. This could be stored locally or in a DB.

    3. The State information. The OpenMM simulation should output a state file
    by using the StateDataReporter. This file should be parsed by OpenMMStateReport
    (below) and saved as an attribute here.

    4. A record of the calculation settings. Each openmm_job should somehow save a
    record (how?) of what happened during the simulation. This should include things like
    the temperature, the ensemble, the number of steps run, and any other settings that
    effect the simulation. Perhaps this should also include the calculation time.
    """


class OpenMMStateReport(BaseModel):
    """OpenMM state."""

    def from_file(self, filename: str):
        """Initialize the state report from an Openmm state output file."""
        pass
