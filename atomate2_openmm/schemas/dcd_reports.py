from pydantic import BaseModel, Field


class DCDReports(BaseModel):
    report_interval: int = Field(None, description="Report interval")
    enforce_periodic_box: bool = Field(None, description="Wrap particles or not")
    blob: bytes = Field(None, description="DCD trajectory blob")
    file_path: str = Field(None, description="file path to document")
    url: str = Field(None, description="URL to locate document")
    host: str = Field(None, description="host name")

    @classmethod
    def from_dcd_file(cls, dcd_file_path, report_interval, enforce_periodic_box):
        with open(dcd_file_path, "rb") as file:
            blob = file.read()
            dcd_reports = DCDReports(
                report_interval=report_interval,
                enforce_periodic_box=enforce_periodic_box,
                blob=blob
            )
            return dcd_reports
