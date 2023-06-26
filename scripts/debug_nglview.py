import nglview as nv
import MDAnalysis as mda


universe = mda.Universe(
    topology="/Users/xperrylinn/Desktop/demo_data/nvt_simple_flow/topology.PDB",
    coordinates="/Users/xperrylinn/Desktop/demo_data/nvt_simple_flow/trajectory.DCD",
)

mda_view = nv.show_mdanalysis(universe.atoms)
mda_view.add_line()
