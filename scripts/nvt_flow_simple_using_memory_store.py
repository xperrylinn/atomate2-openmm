from pymatgen.io.openmm.sets import OpenMMSet
from maggma.stores import MemoryStore
from jobflow import run_locally
from jobflow import JobStore
from jobflow import Flow
from atomate2_openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
from atomate2_openmm.jobs.nvt_maker import NVTMaker

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

# Create JobStore
doc_store = MemoryStore()
trajectory_store = MemoryStore()
job_store = JobStore(docs_store=doc_store, additional_stores={"trajectory_store": trajectory_store})

# Draw the graph before running as it depends on output references
flow.draw_graph().show()

# Run the Production Flow
responses = run_locally(flow=flow, store=job_store, ensure_success=True)

nvt_traj_blob_uuid = next(doc_store.query(criteria={"uuid": flow.jobs[-1].uuid}))["output"]["trajectories"]["blob_uuid"]
dcd_report = next(trajectory_store.query(criteria={"blob_uuid": nvt_traj_blob_uuid}))
assert dcd_report["@class"], "DCDReports"
