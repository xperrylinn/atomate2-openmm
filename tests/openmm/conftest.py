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


@pytest.fixture(scope="function")
def test_output_dir():
    from tempfile import TemporaryDirectory

    # Create a temporary directory to read/write output files
    output_dir = TemporaryDirectory()

    return output_dir.name
