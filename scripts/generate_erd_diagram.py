from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument
import erdantic as erd


diagram = erd.create(OpenMMTaskDocument)

diagram.draw("../schemas_erd.png")
