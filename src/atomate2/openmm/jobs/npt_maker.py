from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from src.atomate2.openmm.schemas.task_details import TaskDetails
from openmm.openmm import MonteCarloBarostat
from dataclasses import dataclass
from openmm.unit import kelvin, atmosphere


@dataclass
class NPTMaker(BaseOpenMMMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    temperature: float = 298
    pressure: float = 1
    frequency: int = 10

    def _run_openmm(self, sim):
        # Add barostat to system
        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        barostat_force_index = system.addForce(
            MonteCarloBarostat(self.pressure * atmosphere, self.temperature * kelvin, 10)
        )

        # Re-init the context after adding thermostat to System
        context.reinitialize(preserveState=True)

        # Run the simulation
        sim.step(self.steps)

        # Remove thermostat and update context
        system.removeForce(barostat_force_index)
        context.reinitialize(preserveState=True)

        task_details = TaskDetails(
            task_name=self.name,
            task_kwargs={
                "steps": self.steps,
                "temperature": self.temperature,
                "frequency": self.frequency,
                "pressure": self.pressure,
            },
            platform_kwargs=self.platform_kwargs,
            total_steps=self.steps,
        )

        return task_details
