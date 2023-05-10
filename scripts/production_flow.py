from src.atomate2.openmm.flows.production_maker import ProductionMaker
from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
from src.atomate2.openmm.jobs.maker_npt import NPTMaker
from pymatgen.io.openmm.sets import OpenMMSet
from maggma.stores import MongoURIStore
from jobflow import run_locally
from jobflow import JobStore
import os


# Define file path to input set of files
input_file_path = "../tests/test_data/alchemy_input_set"
input_set = OpenMMSet.from_directory(
    directory=input_file_path,
    topology_file="topology_pdb",
    system_file="system_xml",
    integrator_file="integrator_xml",
    state_file="state_xml",
    contents_file="contents_json",
)

# Setup a Production Flow
production_maker = ProductionMaker(
    energy_maker=EnergyMinimizationMaker(),
    npt_maker=NPTMaker(steps=100),
    nvt_maker=NVTMaker(steps=100),
)
production_flow = production_maker.make(input_set=input_set)

# Collect database credentials
username, password = os.environ.get("ATLAS_USERNAME"), os.environ.get("ATLAS_PASSWORD")
if username is None or password is None:
    raise ValueError(f"Environment variables ATLAS_USERNAME and ATLAS_PASSWORD must be set.\n"
                     f"ATLAS_USERNAME - Atlas MongoDB username\n"
                     f"ATLAS_PASSWORD - database password\n")

# Setup store using Atlas MongoDB
uri = f"mongodb+srv://{username}:{password}@atomate2-openmm.vlzvqsg.mongodb.net/?retryWrites=true&w=majority"
mongodb_atlas_store = JobStore(
    docs_store=MongoURIStore(
        uri=uri,
        collection_name="Project 0",
        database="atomate2-openmm"
    )
)

# Run the Production Flow
response = run_locally(flow=production_flow, store=mongodb_atlas_store, ensure_success=True)

production_flow.draw_graph().show()