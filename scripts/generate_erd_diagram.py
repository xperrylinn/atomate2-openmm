from atomate2.openmm import OpenMMTaskDocument
import erdantic as erd


diagram = erd.create(OpenMMTaskDocument)

diagram.draw("../schemas_erd.png")
