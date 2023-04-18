from src.atomate2.openmm.flows.core import (
    ProductionMaker,
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

from openmm.app import DCDReporter

from tempfile import (
    NamedTemporaryFile,
    TemporaryDirectory,
)


def test_production_maker_from_input_mol_spec():

    # Define InputMoleculeSpec and other test data
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
    input_mol_dicts = [
        ethanol_molecs,
        water_molecs,
    ]
    density = 0.2

    # Create temp files
    energy_min_tmp_file = NamedTemporaryFile()

    # Define jobs
    input_maker = OpenMMSetFromInputMoleculeSpec()
    energy_maker = EnergyMinimizationMaker()
    npt_maker = NPTMaker(steps=1000)
    anneal_maker = AnnealMaker(steps=[100, 100, 100], temperatures=[500, 1000, 500], temp_steps=[100, 100, 100])
    nvt_maker = NVTMaker(steps=1000)

    production_maker = ProductionMaker(
        name="test_production_maker_from_input_mol_spec",
        input_maker=input_maker,
        energy_maker=energy_maker,
        npt_maker=npt_maker,
        anneal_maker=anneal_maker,
        nvt_maker=nvt_maker,
    )

    # Create a temporary directory for simulation reporters
    temp_dir = TemporaryDirectory()

    flow = production_maker.make(
        output_dir=temp_dir.name,
        input_mol_dicts=input_mol_dicts,
        density=density,
        # box=box,
        default_charge_method="mmff94",
    )

    # flow.draw_graph().show()

    docs_store = MemoryStore()

    store = JobStore(docs_store)

    responses = run_locally(flow, store=store)

    print(responses)


def test_production_maker_from_directory():

    # Define jobs
    input_maker = OpenMMSetFromDirectory()
    energy_maker = EnergyMinimizationMaker()
    npt_maker = NPTMaker(steps=1000)
    anneal_maker = AnnealMaker(
        steps=[100, 100, 100],
        temperatures=[500, 1000, 500],
        temp_steps=[100, 100, 100],
    )
    nvt_maker = NVTMaker(steps=1000)

    # Assumes test are run from repo root
    input_dir = "./tests/data/input_sets/default_input_set_no_dot"

    # Create a temporary directory for simulation reporters
    temp_dir = TemporaryDirectory()

    production_maker = ProductionMaker(
        name="test_production_maker_from_directory",
        input_maker=input_maker,
        energy_maker=energy_maker,
        npt_maker=npt_maker,
        anneal_maker=anneal_maker,
        nvt_maker=nvt_maker,
    )

    flow = production_maker.make(output_dir=temp_dir.name, input_dir=input_dir)

    # flow.draw_graph().show()

    docs_store = MemoryStore()

    store = JobStore(docs_store)

    responses = run_locally(flow, store=store)

    print(responses)
    print("DONE!?!?")


if __name__ == "__main__":
    test_production_maker_from_input_mol_spec()
    test_production_maker_from_directory()
