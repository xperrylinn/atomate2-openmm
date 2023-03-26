from src.atomate2.openmm.flows.core import (
    ProductionMaker,
    ProductionMaker2,
)
from src.atomate2.openmm.jobs.core import (
    OpenMMSetFromDirectory,
    OpenMMSetFromInputMoleculeSpec,
    EnergyMinimizationMaker,
    NPTMaker,
    AnnealMaker,
    NVTMaker,
)

from pymatgen.io.openmm.schema import InputMoleculeSpec

from maggma.stores import MemoryStore

from jobflow.managers.local import run_locally
from jobflow import JobStore


input_molecules = list()
num_water_molecules, num_ethanol_molecules = 400, 20
ethanol_molecs = InputMoleculeSpec(
    name="ethanol",
    smile="CCO",
    count=num_ethanol_molecules,
)
water_molecs = InputMoleculeSpec(
    name="water",
    smile="O",
    count=num_ethanol_molecules,
)
input_molecules.append(ethanol_molecs)
input_molecules.append(water_molecs)

input_mol_dicts = [
    ethanol_molecs,
    water_molecs,
]
density = 1.0
prev_dir = "./data"

# input_maker = OpenMMSetFromInputMoleculeSpec()
input_maker = OpenMMSetFromDirectory()
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

# production_maker = ProductionMaker2(
#     name="my_production_maker",
#     input_maker=input_maker,
#     energy_maker=energy_maker,
#     npt_maker=npt_maker,
#     anneal_maker=anneal_maker,
#     nvt_maker=nvt_maker,
# )

# flow = production_maker.make(
#     input_mol_dicts=input_mol_dicts,
#     density=density,
#     prev_dir=prev_dir,
# )

input_dir = "/Users/xperrylinn/Documents/Academics/MSSE/spring_2023/Capstone/atomate2-openmm/examples/data/input_sets/default_input_set_no_dot"

flow = production_maker.make(input_dir=input_dir)

# flow.draw_graph().show()

docs_store = MemoryStore()

store = JobStore(docs_store)

responses = run_locally(flow, store=store)
