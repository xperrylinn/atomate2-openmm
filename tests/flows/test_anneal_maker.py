def test_anneal_maker(alchemy_input_set, job_store):
    from atomate2_openmm.flows.anneal_maker import AnnealMaker
    from atomate2_openmm.jobs.temp_change_maker import TempChangeMaker
    from atomate2_openmm.jobs.nvt_maker import NVTMaker
    from atomate2_openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from jobflow import run_locally

    anneal_maker = AnnealMaker(
        raise_temp_maker=TempChangeMaker(steps=100, temp_steps=10, final_temp=310),
        nvt_maker=NVTMaker(steps=100, temperature=310),
        lower_temp_maker=TempChangeMaker(steps=100, temp_steps=10),
    )

    anneal_flow = anneal_maker.make(input_set=alchemy_input_set)

    responses = run_locally(flow=anneal_flow, ensure_success=True)

    for job_response in responses.values():
        assert isinstance(job_response[1].output["doc_store"], OpenMMTaskDocument)
