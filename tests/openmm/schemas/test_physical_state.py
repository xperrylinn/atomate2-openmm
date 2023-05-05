def test_phyical_state(test_alchemy_input_set):
    from src.atomate2.openmm.schemas.openmm_task_document import PhysicalState

    physical_state = PhysicalState.from_input_set(test_alchemy_input_set)