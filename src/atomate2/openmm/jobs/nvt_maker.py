from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from src.atomate2.openmm.schemas.task_details import TaskDetails
from openmm.openmm import AndersenThermostat
from dataclasses import dataclass
from openmm.unit import kelvin, pico, second


@dataclass
class NVTMaker(BaseOpenMMMaker):
    steps: int = 1000000
    name: str = "nvt simulation"
    temperature: float = 298

    def _run_openmm(self, sim):
        integrator = sim.context.getIntegrator()
        integrator.setTemperature(self.temperature * kelvin)

        # Run the simulation
        sim.step(self.steps)

        return TaskDetails.from_maker(self)
