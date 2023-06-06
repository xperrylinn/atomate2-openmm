from pymatgen.io.openmm.generators import OpenMMSolutionGen
openmm_sol_gen = OpenMMSolutionGen()
openmm_set = openmm_sol_gen.get_input_set(
    input_mol_dicts=[
        {"smile": "[Li+].[F-]P(F)(F)(F)(F)F", "count": 500},
        {"smile": "CC(=O)OCC", "count": 5000},
        {"smile": "CC(=O)OC(=O)C", "count": 5000},
        {"smile": "CCOC(=O)CC", "count": 200},
        {"smile": "C1C(=O)OC=C1", "count": 100},
        {"smile": "C1C(=O)OCCF1", "count": 100},
    ],
    density=1.0,
)