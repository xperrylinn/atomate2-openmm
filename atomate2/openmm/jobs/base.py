
from typing import Callable
from jobflow import job, Maker
from dataclasses import dataclass
from pathlib import Path
from pymatgen.io.openmm.sets import OpenMMSet
from atomate2.openmm.schemas.task import OpenMMTaskDocument

def openmm_job(method: Callable):
    """

    Parameters
    ----------
    method : callable
        A BaseVaspMaker.make method. This should not be specified directly and is
        implied by the decorator.

    Returns
    -------
    callable
        A decorated version of the make function that will generate OpenMM jobs.
    """
    return job(method, output_schema=OpenMMTaskDocument)


@dataclass
class BaseOpenmmMaker(Maker):
    @openmm_job
    def make(self, input_set: OpenMMSet, prev_dir: str | Path | None = None):
        """
        Run an Openmm calculation.

        Parameters
        ----------
        input_set : OpenMMSet
            A pymatgen structure object.
        prev_dir : str or Path or None
            A previous OpenMM calculation directory to copy output files from.
        """

        # this should mirror the structure of BaseVaspMaker.make
        # it should include everything needed to setup and teardown and openmm
        # simulation but it should offload the details to its children
        # (e.g. NPTMaker, NVTMaker, etc.)

        return
