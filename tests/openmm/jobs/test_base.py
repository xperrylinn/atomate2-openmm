def test_base(test_alchemy_input_set):
    from src.atomate2.openmm.jobs.base import BaseOpenmmMaker
    from tempfile import TemporaryDirectory
    import jobflow

    # Setup
    jstore = jobflow.SETTINGS.JOB_STORE

    base_job_maker = BaseOpenmmMaker()

    base_job = base_job_maker.make(input_set=test_alchemy_input_set)

    with TemporaryDirectory() as temp_dir:

        # Validate _setup_base_openmm_task
        sim = base_job_maker._setup_base_openmm_task(
            input_set=test_alchemy_input_set,
            output_dir=temp_dir,
        )

        # Validate _setup_base_openmm_task
        base_job_maker._close_base_openmm_task(
            input_set=test_alchemy_input_set,
            context=sim.context,
            output_dir=temp_dir,
        )

    # Validate raising of NotImplementedError
    try:
        base_job.run(store=jstore)
        assert False
    except NotImplementedError:
        assert True
