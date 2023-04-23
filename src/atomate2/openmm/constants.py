from enum import Enum


class OpenMMConstants(Enum):
    TOPOLOGY_PDD_FILE_NAME: str = "topology_pdb"
    STATE_XML_FILE_NAME: str = "state_xml"
    SYSTEM_XML_FILE_NAME: str = "system_xml"
    INTEGRATOR_XML_FILE_NAME: str = "integrator_xml"
