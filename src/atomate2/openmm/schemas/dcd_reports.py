from pydantic import BaseModel, Field


class DCDReports(BaseModel):
    report_interval: int = Field(None, description="Report interval")
    enforce_periodic_box: bool = Field(None, description="Wrap particles or not")
    blob: bytes = Field(None, description="DCD trajectory blob")

    @classmethod
    def from_dcd_file(cls, dcd_file):
        # TODO: will somehow need to interface with the additional store?
        return