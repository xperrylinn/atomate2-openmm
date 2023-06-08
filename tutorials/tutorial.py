from pymatgen.io.openmm.schema import InputMoleculeSpec, Geometry
from pymatgen.io.openmm.generators import OpenMMSolutionGen


openmm_set_gen = OpenMMSolutionGen(default_charge_method="mmff94")

partial_charge_scalers = [0.7, 0.8, 0.9, 1.0]

for charge_scaler in partial_charge_scalers:
    ec_mols = InputMoleculeSpec(
        smile="C1COC(=O)O1",
        count=200,
        name="EC",
    )
    emc_mols = InputMoleculeSpec(
        smile="CCOC(=O)OC",
        count=400,
        name="EMC"
    )
    dec_mols = InputMoleculeSpec(
        smile="O=C1OC[C@H](F)O1",
        count=50,
        name="FEC"
    )
    li_mols = InputMoleculeSpec(
        smile="[Li+]",
        count=60,
        name="Li",
        charge_scaling=charge_scaler,
    )
    pf6_mols = InputMoleculeSpec(
        smile="F[P-](F)(F)(F)(F)F",
        count=60,
        name="PF6",
        charge_scaling=charge_scaler,
        geometries=["./data/PF6.xyz"],
        partial_charges=[1.34, -0.39, -0.39, -0.39, -0.39, -0.39, -0.39],
    )

    input_molecules = [ec_mols, emc_mols, dec_mols, li_mols, pf6_mols]

    openmm_set_gen.get_input_set(
        input_mol_dicts=input_molecules,
        density=1.154,
    )

