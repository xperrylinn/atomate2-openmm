from atomate2_openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from atomate2_openmm.schemas.task_details import TaskDetails
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class EnergyMinimizationMaker(BaseOpenMMMaker):
    name: str = "energy minimization"
    # TODO: add default kwargs for Simulation.minimizeEnergy?
    # tolerance
    # maxIterations : int
    dcd_reporter_interval: Optional[int] = field(default=0)
    state_reporter_interval: Optional[int] = field(default=0)

    def _run_openmm(self, sim):

        # Minimize the energy
        sim.minimizeEnergy()
        return TaskDetails.from_maker(self)
