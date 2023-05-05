def test_state_reports(test_state_file):
    from src.atomate2.openmm.schemas.openmm_task_document import StateReports
    state_reporter = StateReports.from_state_file(test_state_file)
