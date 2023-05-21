

with open("data/dcd_test.dcd", "rb") as f:
    data = f.read()


with open("data/dcd_test2.dcd", "wb") as f:
    f.write(data)

import MDAnalysis as mda

u = mda.Universe("data/topology.pdb", "data/dcd_test2.dcd")
print("hi")
