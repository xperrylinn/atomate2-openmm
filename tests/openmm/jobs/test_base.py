def test_base(test_alchemy_input_set):
    from src.atomate2.openmm.jobs.base import BaseOpenmmMaker
    from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from openmm.app.simulation import Simulation
    from tempfile import TemporaryDirectory
    import jobflow

    # Setup
    jstore = jobflow.SETTINGS.JOB_STORE

    base_job_maker = BaseOpenmmMaker(
        state_reporter_interval=0,
        dcd_reporter_interval=0,
    )

    base_job = base_job_maker.make(input_set=test_alchemy_input_set)

    with TemporaryDirectory() as temp_dir:

        # Validate _setup_base_openmm_task
        sim = base_job_maker._setup_base_openmm_task(
            input_set=test_alchemy_input_set,
            output_dir=temp_dir,
        )
        assert isinstance(sim, Simulation)

        # Validate _setup_base_openmm_task
        task_doc = base_job_maker._close_base_openmm_task(
            input_set=test_alchemy_input_set,
            context=sim.context,
            output_dir=temp_dir,
        )
        assert isinstance(task_doc, OpenMMTaskDocument)

    # Validate raising of NotImplementedError
    try:
        base_job.run(store=jstore)
        assert False
    except NotImplementedError:
        assert True
