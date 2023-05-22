def test_nvt_maker(alchemy_input_set, job_store):
    from atomate2_openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from atomate2_openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from atomate2_openmm.jobs.nvt_maker import NVTMaker
    from jobflow import Flow, run_locally

    energy_minimization_job_maker = EnergyMinimizationMaker()
    nvt_job_maker = NVTMaker(
        steps=1000,
        dcd_reporter_interval=10,
        state_reporter_interval=10,
    )
    energy_minimization_job = energy_minimization_job_maker.make(input_set=alchemy_input_set)
    nvt_job = nvt_job_maker.make(input_set=energy_minimization_job.output["doc_store"].calculation_output.input_set)
    flow = Flow(
        jobs=[
            energy_minimization_job,
            nvt_job,
        ]
    )

    response = run_locally(flow=flow, store=job_store, ensure_success=True)
    output = response[nvt_job.uuid][1].output
    assert isinstance(output["doc_store"], OpenMMTaskDocument)

    with job_store.docs_store as s:
        doc = s.query_one({"uuid": nvt_job.uuid})
        assert doc["output"]["doc_store"]["@class"] == "OpenMMTaskDocument"
