from pydantic import BaseModel, Field


class TaskDetails(BaseModel):
    task_name: str = Field(None, description="Task name")
    task_args: dict = Field(None, description="Task args")
    task_kwargs: dict = Field(None, description="Task kwargs")
    platform: str = Field(None, description="Platform")
    platform_details: dict = Field(None, description="Platform details")
    total_steps: int = Field(None, description="Total steps")
