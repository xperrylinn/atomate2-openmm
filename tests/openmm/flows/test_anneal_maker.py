from src.atomate2.openmm.flows.anneal_maker import AnnealMaker
from src.atomate2.openmm.jobs.temperature_maker import TempChangeMaker


def test_anneal_maker(alchemy_input_set, job_store):
    from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
    from src.atomate2.openmm.jobs.npt_maker import NPTMaker
    from src.atomate2.openmm.flows.production_maker import ProductionMaker
    from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
    from jobflow import run_locally

    anneal_maker = AnnealMaker(
        raise_temp_maker=TempChangeMaker(steps=100, temp_steps=10, final_temp=310),
        npt_maker=NPTMaker(steps=100),
        lower_temp_maker=TempChangeMaker(steps=100, temp_steps=10),
    )

    anneal_flow = anneal_maker.make(input_set=alchemy_input_set)

    responses = run_locally(flow=anneal_flow)

    for job_response in responses.values():
        assert isinstance(job_response[1].output, OpenMMTaskDocument)
