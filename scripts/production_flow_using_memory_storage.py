from atomate2_openmm.flows.production_maker import ProductionMaker
from atomate2_openmm.flows.anneal_maker import AnnealMaker
from atomate2_openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
from atomate2_openmm.jobs.nvt_maker import NVTMaker
from atomate2_openmm.jobs.npt_maker import NPTMaker
from atomate2_openmm.jobs.temp_change_maker import TempChangeMaker
from pymatgen.io.openmm.sets import OpenMMSet
from maggma.stores import MongoURIStore
from maggma.stores.aws import S3Store
from maggma.stores import MemoryStore
from jobflow import run_locally
from jobflow import JobStore
import os


# Define file path to input set of files
input_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../tests/test_data/alchemy_input_set")

# Create the initial OpenMMSet
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
    npt_maker=NPTMaker(
        steps=100,
        state_reporter_interval=10,
        dcd_reporter_interval=10,
    ),
    anneal_maker=AnnealMaker(
        raise_temp_maker=TempChangeMaker(
            steps=1000,
            temp_steps=10,
            final_temp=700,
            state_reporter_interval=0,
            dcd_reporter_interval=0,
        ),
        nvt_maker=NVTMaker(
            steps=100,
            state_reporter_interval=0,
            dcd_reporter_interval=0,
            temperature=700,
        ),
        lower_temp_maker=TempChangeMaker(
            steps=1000,
            temp_steps=100,
            final_temp=298,
            state_reporter_interval=0,
            dcd_reporter_interval=0,
        ),
    ),
    nvt_maker=NVTMaker(
        steps=100,
        state_reporter_interval=10,
        dcd_reporter_interval=10,
    ),
)
production_flow = production_maker.make(input_set=input_set)

# Setup memory store for document store
doc_store = MemoryStore()

# Setup memory store for reporter store
trajectory_store = MemoryStore()

# Create JobStore
job_store = JobStore(
    docs_store=doc_store,
    additional_stores={"trajectory_store": trajectory_store},
)

# Draw the graph before running as it depends on output references
production_flow.draw_graph().show()

# Run the Production Flow
response = run_locally(flow=production_flow, store=job_store, ensure_success=True)
