import time
import pathlib
import sys
import os

from openmm.app import DCDReporter, StateDataReporter, PDBReporter
from pymatgen.io.openmm.sets import OpenMMSet
from openmm import Platform
from pymatgen.io.openmm.simulations import anneal, equilibrate_pressure


"""
argument 1: the input directory
argument 2: the output directory
"""


wd = pathlib.Path.cwd()
input_dir = pathlib.Path(sys.argv[1])
output_dir = pathlib.Path(sys.argv[2])

target_paths = [path for path in input_dir.iterdir()]
assert len(target_paths) == 1, "there must be exactly 1 job to run"

target_path = target_paths[0]
input_set = OpenMMSet.from_directory(target_path)
print("OpenMMSet created")

platform = Platform.getPlatformByName("CPU")
sim = input_set.get_simulation()
print("OpenMM simulation created")

output_dir = output_dir / target_path.stem
output_dir.mkdir(parents=True, exist_ok=True)
pdb_reporter = PDBReporter(str(output_dir / "topology_pdb"), 1)
pdb_reporter.report(sim, sim.context.getState(getPositions=True))
sim.reporters.append(
    StateDataReporter(
        str(output_dir / "state_txt"),
        10,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True,
    )
)
sim.reporters.append(DCDReporter(str(output_dir / "trajectory_equil_dcd"), 10))
print("StateDataReporter addedto simulation reporters")

start = time.time()

sim.minimizeEnergy()
print("Simulation particle coordinates energy minimized")

equilibrate_pressure(sim, 100)
print("Simulation pressure equilibrated")

anneal(sim, 400, [100, 100, 100])
print("Simulation annealed")

sim.reporters.pop()
sim.reporters.append(DCDReporter(str(output_dir / "trajectory_dcd"), 10))
print("Simulation DCDReporter added")

sim.step(5000)
print("Simulation step")

print("total runtime: ", time.time() - start)
