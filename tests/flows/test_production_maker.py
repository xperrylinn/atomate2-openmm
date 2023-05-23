def test_production_maker(alchemy_input_set, job_store):
    from atomate2_openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from atomate2_openmm.jobs.nvt_maker import NVTMaker
    from atomate2_openmm.jobs.npt_maker import NPTMaker
    from atomate2_openmm.flows.production_maker import ProductionMaker
    from atomate2_openmm.flows.anneal_maker import AnnealMaker
    from atomate2_openmm.schemas.openmm_task_document import OpenMMTaskDocument
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

    production_flow = production_maker.make(input_set=alchemy_input_set)

    responses = run_locally(flow=production_flow, ensure_success=True)

    for job_response in responses.values():
        assert isinstance(job_response[1].output["doc_store"], OpenMMTaskDocument)
