import MDAnalysis as mda
import io

topology_file_path = "/Users/xperrylinn/Desktop/demo_data/3_nvt_simulation/topology.pdb"
trajectory_file_path = "/Users/xperrylinn/Desktop/demo_data/3_nvt_simulation/trajectory.dcd"

with open(trajectory_file_path, 'rb') as file:
    trajectory_blob = file.read()

str_stream = io.StringIO(topology_file_path)
bytes_stream = io.BytesIO(trajectory_blob)

universe = mda.Universe(
    topology=topology_file_path,
    topology_format="PDB",
    coordinates=bytes_stream,
    format="DCD",
)