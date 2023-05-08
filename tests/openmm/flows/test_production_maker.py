def test_production_maker(alchemy_input_set, job_store):
    from src.atomate2.openmm.jobs.energy_minimization_maker import EnergyMinimizationMaker
    from src.atomate2.openmm.jobs.nvt_maker import NVTMaker
    from src.atomate2.openmm.jobs.maker_npt import NPTMaker
    from src.atomate2.openmm.flows.production_maker import ProductionMaker
    from jobflow import run_locally

    production_maker = ProductionMaker(
        energy_maker=EnergyMinimizationMaker(),
        npt_maker=NPTMaker(steps=100),
        nvt_maker=NVTMaker(steps=100),
    )

    production_flow = production_maker.make(input_set=alchemy_input_set)

    responses = run_locally(flow=production_flow, ensure_success=True)
