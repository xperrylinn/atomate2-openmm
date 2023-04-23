def test_energy_minimization_maker(test_input_dir, test_output_dir):
    import jobflow
    from jobflow import run_locally
    from src.atomate2.openmm.jobs.core import EnergyMinimizationMaker
    from src.atomate2.openmm.constants import OpenMMConstants
    from pymatgen.io.openmm.sets import OpenMMSet

    jstore = jobflow.SETTINGS.JOB_STORE

    # Create OpenMM set from directory for job input
    input_set = OpenMMSet.from_directory(
        directory=test_input_dir,
        topology_file=OpenMMConstants.TOPOLOGY_PDD_FILE_NAME.value,
        state_file=OpenMMConstants.STATE_XML_FILE_NAME.value,
        system_file=OpenMMConstants.SYSTEM_XML_FILE_NAME.value,
        integrator_file=OpenMMConstants.INTEGRATOR_XML_FILE_NAME.value,
    )

    # Instantiate a Maker and make a Job
    job = EnergyMinimizationMaker().make(input_set=input_set, output_dir=test_output_dir)

    # Run the the job
    responses = run_locally(flow=job, ensure_success=True)

    # Assertions on result
    output = responses[job.uuid][1].output
    assert isinstance(output, OpenMMSet)
    # todo: mock or not-mock and assert PE is < initial PE?

    with jstore.additional_stores["data"] as s:
        doc = s.query_one({"job_uuid": job.uuid})
        dd = doc["data"]
        assert dd["@class"] == "Chgcar"


def test_npt_maker():
    pass


def test_nvt_maker():
    pass


def test_anneal_maker():
    pass
