from pymatgen.io.openmm.sets import OpenMMSet
from src.atomate2.openmm.constants import OpenMMConstants
import pytest


@pytest.fixture(scope="session")
def platform_for_test():
    return "CPU"


@pytest.fixture(scope="session")
def platform_props_for_test():
    return None


@pytest.fixture(scope="session")
def test_input_dir():
    import os
    from pathlib import Path

    module_dir = Path(__file__).resolve().parent
    test_dir = Path(os.path.join(module_dir, "../test_data/alchemy_input_set"))

    return test_dir.resolve()


@pytest.fixture(scope="session")
def test_input_set(test_input_dir):

    input_set = OpenMMSet.from_directory(
        directory=test_input_dir,
        topology_file=OpenMMConstants.TOPOLOGY_PDD_FILE_NAME.value,
        state_file=OpenMMConstants.STATE_XML_FILE_NAME.value,
        system_file=OpenMMConstants.SYSTEM_XML_FILE_NAME.value,
        integrator_file=OpenMMConstants.INTEGRATOR_XML_FILE_NAME.value,
        contents_file=OpenMMConstants.CONTENTS_JOSN_FILE_NAME.value,
    )
    return input_set


@pytest.fixture
def test_state_file(test_input_dir):
    import os
    from pathlib import Path

    state_file = Path(os.path.join(test_input_dir, "state_xml"))

    return state_file.resolve()


@pytest.fixture(scope="function")
def test_output_dir(tmp_path_factory):
    return tmp_path_factory.mktemp(basename="output")
