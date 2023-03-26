from src.atomate2.openmm.jobs.core import (
    OpenMMSetFromInputMoleculeSpec,
    OpenMMSetFromDirectory,
    EnergyMinimizationMaker,
    NPTMaker,
    AnnealMaker,
    NVTMaker,
)
from typing import List, Union, Optional, Dict
from pathlib import Path

from dataclasses import dataclass

from pydantic import Field
from jobflow import Maker, Flow, job

from pymatgen.io.openmm.schema import InputMoleculeSpec


@dataclass
class ProductionMaker(Maker):
    """


    """
    name: str = "production"

    input_maker: OpenMMSetFromDirectory = Field(default_factory=OpenMMSetFromDirectory)
    energy_maker: EnergyMinimizationMaker = Field(default_factory=EnergyMinimizationMaker)
    npt_maker: NPTMaker = Field(default_factory=NPTMaker)
    anneal_maker: AnnealMaker = Field(default_factory=AnnealMaker)
    nvt_maker: NVTMaker = Field(default_factory=NVTMaker)

    def make(self, input_dir: Union[str, Path]):

        input_job = self.input_maker.make(input_dir=input_dir)

        energy_job = self.energy_maker.make(input_set=input_job.output)

        pressure_job = self.npt_maker.make(input_set=energy_job.output.input_set)

        anneal_job = self.anneal_maker.make(input_set=pressure_job.output.input_set)

        production_job = self.nvt_maker.make(input_set=anneal_job.output.input_set)

        my_flow = Flow(
            [input_job, energy_job, pressure_job, anneal_job, production_job],
            output={"log": production_job},
        )

        return my_flow


@dataclass
class ProductionMaker2(Maker):
    """


    """
    name: str = "production2"

    input_maker: OpenMMSetFromInputMoleculeSpec = Field(default_factory=OpenMMSetFromInputMoleculeSpec)
    energy_maker: EnergyMinimizationMaker = Field(default_factory=EnergyMinimizationMaker)
    npt_maker: NPTMaker = Field(default_factory=NPTMaker)
    anneal_maker: AnnealMaker = Field(default_factory=AnnealMaker)
    nvt_maker: NVTMaker = Field(default_factory=NVTMaker)

    def make(self, input_dir):
        input_job = self.input_maker.make(input_dir=input_dir)

        energy_job = self.energy_maker.make(input_set=input_job.output.input_set)

        pressure_job = self.npt_maker.make(input_set=energy_job.output.input_set)

        anneal_job = self.anneal_maker.make(input_set=pressure_job.output.input_set)

        production_job = self.nvt_maker.make(input_set=anneal_job.output.input_set)

        my_flow = Flow(
            [input_job, energy_job, pressure_job, anneal_job, production_job],
            output={"log": production_job},
        )

        return my_flow
