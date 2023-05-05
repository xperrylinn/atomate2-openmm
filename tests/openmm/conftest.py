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
def test_data_dir():
    import os
    from pathlib import Path

    module_dir = Path(__file__).resolve().parent.parent
    test_data_dir = Path(os.path.join(module_dir, "./test_data")).resolve()

    return test_data_dir


@pytest.fixture(scope="session")
def test_alchemy_input_set_dir(test_data_dir):
    import os
    from pathlib import Path

    alchemy_input_set_dir = Path(os.path.join(test_data_dir, "./alchemy_input_set")).resolve()

    return alchemy_input_set_dir


@pytest.fixture(scope="session")
def test_alchemy_input_set(test_alchemy_input_set_dir):

    input_set = OpenMMSet.from_directory(
        directory=test_alchemy_input_set_dir,
        topology_file=OpenMMConstants.TOPOLOGY_PDD_FILE_NAME.value,
        state_file=OpenMMConstants.STATE_XML_FILE_NAME.value,
        system_file=OpenMMConstants.SYSTEM_XML_FILE_NAME.value,
        integrator_file=OpenMMConstants.INTEGRATOR_XML_FILE_NAME.value,
        contents_file=OpenMMConstants.CONTENTS_JOSN_FILE_NAME.value,
    )
    return input_set


@pytest.fixture
def test_state_report_file(test_data_dir):
    import os
    from pathlib import Path

    state_file = Path(os.path.join(test_data_dir, "./reporters/state.txt"))

    return state_file.resolve()


@pytest.fixture(scope="function")
def test_output_dir(tmp_path_factory):
    return tmp_path_factory.mktemp(basename="output")
