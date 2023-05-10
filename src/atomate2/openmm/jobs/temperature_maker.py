from src.atomate2.openmm.jobs.base_openmm_maker import BaseOpenMMMaker
from src.atomate2.openmm.schemas.task_details import TaskDetails
from openmm.openmm import MonteCarloBarostat
from dataclasses import dataclass
from openmm.unit import kelvin, atmosphere
import numpy as np


@dataclass
class TempChangeMaker(BaseOpenMMMaker):
    steps: int = 1000000
    name: str = "npt simulation"
    final_temp: float = 298
    temp_steps: int = 100

    def _run_openmm(self, sim):
        # Add barostat to system
        integrator = sim.context.getIntegrator()
        start_temp = integrator.getTemperature()

        # Heating temperature
        delta_t = abs(self.final_temp * kelvin - start_temp)
        temp_step_size = delta_t / self.temp_steps
        for temp in np.arange(
                start_temp + temp_step_size,
                self.final_temp * kelvin + temp_step_size,
                temp_step_size,
        ):
            integrator.setTemperature(temp * kelvin)
            sim.step(self.steps // self.temp_steps)

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
