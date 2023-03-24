from src.atomate2.openmm.flows.core import ProductionMaker
from src.atomate2.openmm.jobs.core import (
    InputMaker,
    EnergyMinimizationMaker,
    NPTMaker,
    AnnealMaker,
    NVTMaker,
)
from pymatgen.io.openmm.schema import InputMoleculeSpec


input_molecules = list()
num_water_molecules, num_ethanol_molecules = 400, 20
ethanol_molecs = InputMoleculeSpec(name="ethanol", smile="O", count=num_ethanol_molecules)
water_molecs = InputMoleculeSpec(name="ethanol", smile="O", count=num_ethanol_molecules)
input_molecules.append(ethanol_molecs)
input_molecules.append(water_molecs)

input_mol_dicts = [
    ethanol_molecs,
    water_molecs,
]
density = 1.0
prev_dir = "./data"

input_maker = InputMaker()
energy_maker = EnergyMinimizationMaker()
npt_maker = NPTMaker()
anneal_maker = AnnealMaker()
nvt_maker = NVTMaker()

production_maker = ProductionMaker(
    name="my_production_maker",
    input_maker=input_maker,
    energy_maker=energy_maker,
    npt_maker=npt_maker,
    anneal_maker=anneal_maker,
    nvt_maker=nvt_maker,
)

flow = production_maker.make(
    input_mol_dicts=input_mol_dicts,
    density=density,
    prev_dir=prev_dir,
)

flow.draw_graph().show()