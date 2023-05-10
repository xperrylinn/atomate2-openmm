from src.atomate2.openmm.flows.anneal_maker import AnnealMaker


def test_production_maker(alchemy_input_set, job_store):
    from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
    from src.atomate2.openmm.jobs.npt_maker import NPTMaker
    from src.atomate2.openmm.flows.production_maker import ProductionMaker
    from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from jobflow import run_locally

    anneal_maker = AnnealMaker.from_temps_and_steps(
        anneal_temp=310,
        final_temp=298,
        steps=100,
        temp_steps=10
    )

    production_maker = ProductionMaker(
        energy_maker=EnergyMinimizationMaker(),
        npt_maker=NPTMaker(steps=100),
        anneal_maker=anneal_maker,
        nvt_maker=NVTMaker(steps=100),
    )

    production_flow = production_maker.make(input_set=alchemy_input_set)

    responses = run_locally(flow=production_flow)

    for job_response in responses.values():
        assert isinstance(job_response[1].output, OpenMMTaskDocument)
