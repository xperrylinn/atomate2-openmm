from pymatgen.io.openmm.sets import OpenMMSet
from maggma.stores import MongoURIStore
from maggma.stores.aws import S3Store
from maggma.stores import MemoryStore
from jobflow import run_locally
from jobflow import JobStore
from jobflow import Flow
from atomate2_openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
from atomate2_openmm.jobs.nvt_maker import NVTMaker
from atomate2_openmm.schemas.dcd_reports import DCDReports
from tempfile import TemporaryDirectory
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

# Create jobs
energy_minimization_job = EnergyMinimizationMaker().make(input_set=input_set)
nvt_job = NVTMaker(
    steps=100,
    state_reporter_interval=10,
    dcd_reporter_interval=10,
    temperature=700,
).make(input_set=energy_minimization_job.output["doc_store"].calculation_output.output_set)

# Setup a Flow
flow = Flow(jobs=[energy_minimization_job, nvt_job],)

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

# Setup S3Store for reporter files
temp_dir = TemporaryDirectory()
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

# Run the Production Flow
responses = run_locally(flow=flow, store=job_store, ensure_success=True)

responses[flow.jobs[-1].uuid][1].output["trajectories"]
responses[flow.jobs[-1].uuid][1].output["doc_store"]

nvt_traj_blob_uuid = next(atlas_mongo_store.query(criteria={"uuid": flow.jobs[-1].uuid}))["output"]["trajectories"]["blob_uuid"]
dcd_report = next(s3_store.query(criteria={"blob_uuid": nvt_traj_blob_uuid}))
assert dcd_report["@class"], "DCDReports"
