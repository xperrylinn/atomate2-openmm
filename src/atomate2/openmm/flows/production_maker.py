from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
from src.atomate2.openmm.jobs.npt_maker import NPTMaker
from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
from src.atomate2.openmm.flows.anneal_maker import AnnealMaker
from pymatgen.io.openmm.sets import OpenMMSet
from typing import Optional, Union
from pathlib import Path
from dataclasses import dataclass, field
from pydantic import Field
from jobflow import Maker, Flow


@dataclass
class ProductionMaker(Maker):
    """
    Class for running
    """
    name: str = "production"
    energy_maker: EnergyMinimizationMaker = field(default_factory=EnergyMinimizationMaker)
    npt_maker: NPTMaker = field(default_factory=NPTMaker)
    anneal_maker: AnnealMaker = field(default_factory=AnnealMaker)
    nvt_maker: NVTMaker = field(default_factory=NVTMaker)

    def make(self, input_set: OpenMMSet, output_dir: Optional[Union[str, Path]] = None):
        """

        Parameters
        ----------
        input_set : OpenMMSet
            OpenMMSet object instance.
        output_dir : Optional[Union[str, Path]]
            File path to directory for writing simulation output files (e.g. trajectory and state files)

        Returns
        -------

        """

        energy_job = self.energy_maker.make(
            input_set=input_set,
            output_dir=output_dir
        )

        pressure_job = self.npt_maker.make(
            input_set=energy_job.output["doc_store"].calculation_output.output_set,
            output_dir=output_dir
        )

        anneal_job = self.anneal_maker.make(
            input_set=pressure_job.output["doc_store"].calculation_output.output_set,
            output_dir=output_dir
        )

        nvt_job = self.nvt_maker.make(
            input_set=anneal_job.output["doc_store"].calculation_output.output_set,
            output_dir=output_dir
        )

        my_flow = Flow(
            [
                energy_job,
                pressure_job,
                anneal_job,
                nvt_job,
            ],
            output={"log": nvt_job, "test": anneal_job},
        )

        return my_flow
