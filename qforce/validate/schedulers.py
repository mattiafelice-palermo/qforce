# A scheduler file should be composed of the following parts
# 1. Scheduler settings
#  - Nodes
#  - number of tasks
#  - ram specification
#  - queue
#  - etc.
# 2. Preparatory phase
#  - module activation
#  - conda activation
#  - handling of scratch dir
#  - etc
#  3. Run command
#  - Can either use default command or user customized commands
#  4. Finalization
#  - handling of scratch dir
#  - cleanup
#  - send signal for any dependent job?
#  - etc

from colt import Colt
from abc import ABC, abstractmethod
import shutil
import os

from pprint import pprint


def get_scheduler(settings):
    """
    Returns the appropriate scheduler based on the specified settings.

    Args:
        settings: The settings object containing the configurations for qforce-validate

    Returns:
        Scheduler: An instance of a scheduler based on the specified options.

    Raises:
        NotImplementedError: If the scheduler specified in the settings is not supported.
    """
    scheduler = settings.general.scheduler.lower()

    if scheduler == "manual":
        return ManualScheduler(settings)
    elif scheduler == "system":
        return SystemScheduler(settings)
    elif scheduler == "pbs":
        pass
    elif scheduler == "slurm":
        pass
    else:
        raise NotImplementedError(f"Scheduler '{scheduler}' is not implemented.")


class SchedulerABC(ABC, Colt):
    _user_input = """
    # 
    n_tasks = 8 :: int

    # In gigabytes - total
    ram_size = 10 :: int
    
    #
    queue_name = :: str, optional

    # scratch directory
    scratch_dir = :: str, optional

    """

    def __init__(self, settings):
        self.settings = settings
        self._jobs = []
        self._n_tasks = settings.scheduler.n_tasks
        self._ram_size = settings.scheduler.ram_size
        self._queue_name = settings.scheduler.queue_name
        self._scratch_dir = settings.scheduler.scratch_dir

    def add(self, job):
        self._jobs.append(job)

    @abstractmethod
    def execute(self): ...

    @property
    def n_tasks(self):
        return self._n_tasks

    @property
    def ram_size(self):
        return self._ram_size

    @property
    def queue_name(self):
        return self._queue_name

    @property
    def scratch_dir(self):
        return self._scratch_dir


class ManualScheduler(SchedulerABC):
    """Just instruct generators and resamplers to write input files and launch scripts to disk."""

    def execute(self):
        for job in self._jobs:
            job.run(dry_run=True)


class SystemScheduler(SchedulerABC):
    """Just instruct generators and resamplers to write input files and launch scripts to disk."""

    def execute(self):
        for job in self._jobs:
            job.run()
