from pydantic import BaseModel, Field


class TaskDetails(BaseModel):
    task_name: str = Field(None, description="Task name")
    task_kwargs: dict = Field(None, description="Task kwargs")
    platform_kwargs: dict = Field(None, description="Platform")
    total_steps: int = Field(None, description="Total steps")
