from jobflow import run_locally

from atomate2.openmm.jobs.openmm_set_maker import OpenMMSetMaker

from maggma.stores import MemoryStore

def test_openmm_set_maker(job_store):
    input_mol_dicts = [
        {"smile": "O", "count": 200},
        {"smile": "CCO", "count": 20},
    ]
    from atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from atomate2.openmm.jobs.nvt_maker import NVTMaker
    from atomate2.openmm.jobs.npt_maker import NPTMaker
    from atomate2.openmm.flows.production_maker import ProductionMaker
    from atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from jobflow import run_locally

    anneal_maker = AnnealMaker.from_temps_and_steps(
        anneal_temp=310,
        final_temp=298,
        steps=100,
        temp_steps=10,
        base_kwargs=dict(
            state_reporter_interval=10,
            dcd_reporter_interval=10,
        )
    )
    energy_maker = EnergyMinimizationMaker()
    npt_maker = NPTMaker(
        steps=100,
        state_reporter_interval=10,
        dcd_reporter_interval=10,
    )
    nvt_maker = NVTMaker(
        steps=100,
        state_reporter_interval=10,
        dcd_reporter_interval=10,
    )

    production_maker = ProductionMaker(
        energy_maker=energy_maker,
        npt_maker=npt_maker,
        anneal_maker=anneal_maker,
        nvt_maker=nvt_maker,
    )

    set_maker = OpenMMSetMaker()

    set_job = set_maker.make(input_mol_dicts=input_mol_dicts, density=1)

    responses = run_locally(set_job, store=job_store)

    print("hi")
