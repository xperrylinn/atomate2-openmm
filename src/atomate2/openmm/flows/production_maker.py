from src.atomate2.openmm.jobs.openmm_set_from_input_mol_spec import (
    OpenMMSetFromInputMoleculeSpec,
    OpenMMSetFromDirectory,
    EnergyMinimizationMaker,
    NPTMaker,
    AnnealMaker,
    NVTMaker,
)
from typing import (
    List,
    Union,
    Optional,
    Dict,
    Union,
)
from pathlib import Path

from dataclasses import dataclass

from pydantic import Field
from jobflow import Maker, Flow, job

from pymatgen.io.openmm.schema import InputMoleculeSpec


@dataclass
class ProductionMaker(Maker):
    """
    Class for running
    """
    name: str = "production"

    input_maker: OpenMMSetFromDirectory = Field(default_factory=OpenMMSetFromDirectory)
    energy_maker: EnergyMinimizationMaker = Field(default_factory=EnergyMinimizationMaker)
    npt_maker: NPTMaker = Field(default_factory=NPTMaker)
    anneal_maker: AnnealMaker = Field(default_factory=AnnealMaker)
    nvt_maker: NVTMaker = Field(default_factory=NVTMaker)

    def make(self, input_mol_spec: InputMoleculeSpec, output_dir: Optional[Union[str, Path]]):
        """

        Parameters
        ----------
        output_dir : Optional[Union[str, Path]]
            File path to directory for writing simulation output files (e.g. trajectory and state files)

        Returns
        -------

        """

        energy_job = self.energy_maker.make(input_set=openmm_set_gen.output, output_dir=output_dir)

        pressure_job = self.npt_maker.make(input_set=energy_job.output, output_dir=output_dir)

        anneal_job = self.anneal_maker.make(input_set=pressure_job.output, output_dir=output_dir)

        production_job = self.nvt_maker.make(input_set=anneal_job.output, output_dir=output_dir)

        my_flow = Flow(
            [
                energy_job,
                pressure_job,
                anneal_job,
                production_job,
            ],
            output={"log": production_job},
        )

        return my_flow
