def test_task_details(platform, platform_properties):
    from atomate2_openmm.schemas.task_details import TaskDetails

    task_details = TaskDetails(
        task_name="my_task",
        task_kwargs={
            "temperature": 100.0,
            "steps": 1000,
            "frequency": 10,
        },
        platform_kwargs={
            "platform": platform,
            "platform_properties": platform_properties,
        },
        total_steps=1000,
    )
