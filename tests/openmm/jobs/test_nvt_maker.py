def test_nvt_maker(alchemy_input_set, job_store):
    from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
    from jobflow import Flow, run_locally

    energy_minimization_job_maker = EnergyMinimizationMaker()
    nvt_job_maker = NVTMaker(
        steps=1000,
        dcd_reporter_interval=10,
        state_reporter_interval=10,
    )
    energy_minimization_job = energy_minimization_job_maker.make(input_set=alchemy_input_set)
    nvt_job = nvt_job_maker.make(input_set=energy_minimization_job.output.calculation_output.output_set)
    flow = Flow(
        jobs=[
            energy_minimization_job,
            nvt_job,
        ]
    )
    run_locally(flow=flow, store=job_store, ensure_success=True)
