from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenmmMaker
from src.atomate2.openmm.schemas.task_details import TaskDetails
from openmm.openmm import AndersenThermostat
from dataclasses import dataclass
from openmm.unit import kelvin, pico, second


@dataclass
class NVTMaker(BaseOpenmmMaker):
    steps: int = 1000000
    name: str = "nvt simulation"
    temperature: float = 298
    frequency: int = 10

    def _run_openmm(self, sim):    # parameters should be NVTMaker attrs
        context = sim.context
        system = context.getSystem()
        assert (
            system.usesPeriodicBoundaryConditions()
        ), "system must use periodic boundary conditions for pressure equilibration."
        thermostat_force_index = system.addForce(
            AndersenThermostat(self.temperature * kelvin, pico * second * self.frequency)
        )

        # Re-init the context afte adding thermostat to System
        context.reinitialize(preserveState=True)

        # Run the simulation
        sim.step(self.steps)

        # Remove thermostat and update context
        system.removeForce(thermostat_force_index)
        context.reinitialize(preserveState=True)

        # TODO: return TaskDetails
        task_details = TaskDetails(
            task_name=self.name,
            task_kwargs={
                "steps"=self.steps,
                "temperature"=self.temperature,
                "frequency"=self.frequency,
            }
        )

        return task_details
