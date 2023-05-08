from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from src.atomate2.openmm.schemas.task_details import TaskDetails
from dataclasses import dataclass


@dataclass
class EnergyMinimizationMaker(BaseOpenMMMaker):
    name: str = "energy minimization"
    # TODO: add default kwargs for Simulation.minimizeEnergy?
    # tolerance
    # maxIterations : int
    state_reporter_interval: int = 0

    def _run_openmm(self, sim):

        # Minimize the energy
        sim.minimizeEnergy()

        task_details = TaskDetails(
            task_name=self.name,
            task_kwargs=None,
            platform_kwargs=self.platform_kwargs,
            total_steps=0,
        )

        return task_details
