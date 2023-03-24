from typing import Optional, Union, Callable
from jobflow import job, Maker
from dataclasses import dataclass
from pathlib import Path
from pymatgen.io.openmm.sets import OpenMMSet
from src.atomate2.openmm.schemas.openmm_task_document import OpenMMTaskDocument


def openmm_job(method: Callable):
    """

    Parameters
    ----------
    method : callable
        A BaseOpenmmMaker.make method. This should not be specified directly and is
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
    def make(self, input_set: OpenMMSet, prev_dir: Optional[Union[str, Path]] = None):
        """
        Run an Openmm calculation.

        Parameters
        ----------
        input_set : OpenMMSet
            pymatgen.io.openmm OpenMMSet object instance.
        prev_dir : str or Path or None
            Previous OpenMM calculation directory to copy output files from.
        """

        # this should mirror the structure of BaseVaspMaker.make
        # it should include everything needed to setup and teardown and openmm
        # simulation but it should offload the details to its children
        # (e.g. NPTMaker, NVTMaker, etc.)

        return
