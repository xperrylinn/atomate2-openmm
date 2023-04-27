
import pytest
import pathlib

from src.atomate2.openmm.schemas.openmm_task_document import StateReports


def test_from_state_file():
    # TODO: not use a garbage hard coded string
    state_file = "/Users/orioncohen/projects/development/atomate2-openmm/tests/test_data/reporters/state.txt"
    with open(state_file, 'r') as file:
        state_str = file.read()
    state_file = pathlib.Path(state_file)
    reporter = StateReports.from_state_file(state_file)
    return
