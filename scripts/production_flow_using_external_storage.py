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

# Collect Atlas database credentials
username, password = os.environ.get("ATLAS_USERNAME"), os.environ.get("ATLAS_PASSWORD")
if username is None or password is None:
    raise ValueError(f"Environment variables ATLAS_USERNAME and ATLAS_PASSWORD must be set.\n"
                     f"ATLAS_USERNAME - Atlas MongoDB username\n"
                     f"ATLAS_PASSWORD - database password\n")

# Setup store using Atlas MongoDB
uri = f"mongodb+srv://{username}:{password}@atomate2-openmm.vlzvqsg.mongodb.net/?retryWrites=true&w=majority"
atlas_mongo_store = MongoURIStore(
    uri=uri,
    collection_name="Project 0",
    database="atomate2-openmm"
)

# Setup S3 store for reporter blobs
index = MemoryStore(collection_name="index", key="blob_uuid")
s3_store = S3Store(
    index=index,
    bucket="atomate2-openmm",
    endpoint_url="https://s3.us-west-1.amazonaws.com",
    s3_profile="atomate2-openmm-dev",
    key="blob_uuid",
    s3_workers=1,
    unpack_data=True,
)

# Create JobStore
job_store = JobStore(
    docs_store=atlas_mongo_store,
    additional_stores={"trajectory_store": s3_store},
)

# Draw the graph before running as it depends on output references
production_flow.draw_graph().show()

# Run the Production Flow
response = run_locally(flow=production_flow, store=job_store, ensure_success=True)