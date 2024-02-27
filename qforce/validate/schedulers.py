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
import textwrap
import subprocess
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
        return PbsScheduler(settings)
    elif scheduler == "slurm":
        return SlurmScheduler(settings)
    else:
        raise NotImplementedError(f"Scheduler '{scheduler}' is not implemented.")


class SchedulerABC(ABC, Colt):
    _user_input = """
    #
    queue_name = :: str, optional

    # scratch directory
    scratch_dir = :: str, optional

    """

    def __init__(self, settings):
        self.settings = settings
        self._pending_jobs = []
        self._completed_jobs = []
        self._queue_name = settings.scheduler.queue_name
        self._scratch_dir = settings.scheduler.scratch_dir

    def add(self, job):
        self._pending_jobs.append(job)

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
        for job in self._pending_jobs:
            job.run(dry_run=True)


class SystemScheduler(SchedulerABC):
    """Just instruct generators and resamplers to write input files and launch scripts to disk."""

    def execute(self):
        for job in self._pending_jobs:
            job.run()


class SlurmScheduler(SchedulerABC):
    def __init__(Self, settings):
        super().__init__(settings)

    def execute(self):
        job_ids = []

        if not all([job.queue for job in self._pending_jobs]):
            # TODO: provide a more specific error pointing to which job is lacking a queue
            raise ValueError(
                "Jobs are missing a work queue."
                "Please make sure that all generators and calculators have a set work queue when using the Slurm scheduler"
            )

        for job in self._pending_jobs:
            job.run(dry_run=True)
            job_dir = job.job_dir

            with open(os.path.join(job_dir, "submit_job.slm"), "w") as slurm_file_handle:
                slurm_file_handle.write(self._slurm_script_content(job))

            job_id = subprocess.run(
                "sbatch submit_job.slm | awk '{print $4}'", cwd=job_dir, shell=True, capture_output=True, text=True
            )
            print(job_id)
            job_ids.append(job_id.stdout.strip("\n"))

        job_string = ":".join([f"{job_id}" for job_id in job_ids])

        print(job_string)
        dummy_job = subprocess.run(f"sbatch --wait --dependency=afterany:{job_string} --wrap='sleep 1'", shell=True)

        # Once completed, move the jobs out of the pending job list
        while self._pending_jobs:
            self._completed_jobs.append(self._pending_jobs.pop())

    def _slurm_script_content(self, job):

        scratch_dir = os.path.join(self.scratch_dir, "$USER")
        scratch_partition_management_pre = f"""
            # Create user-specific scratch directory if it does not exist
            if ! [ -d {scratch_dir} ]; then
            mkdir {scratch_dir}
            fi

            # Create a job-specific scratch directory and set environment variable
            export SCRDIR={scratch_dir}/$SLURM_JOB_ID
            mkdir $SCRDIR

            ## Go to scratch directory and copy everything from working directory
            cd $SCRDIR
            cp -pr $SLURM_SUBMIT_DIR/* . \n
            """

        scratch_partition_management_post = f"""
            # Copy results back to the original working directory
            cp -pr * $SLURM_SUBMIT_DIR/

            # Cleanup scratch directory after successful job completion
            cd $SLURM_SUBMIT_DIR
            rm -rf $SCRDIR
            """

        if self.scratch_dir is None:
            scratch_partition_management_pre = ""
            scratch_partition_management_post = ""

        if job.total_memory is None:
            job_memory_string = ""
        else:
            job_memory_string = f"#SBATCH --mem= {job.total_memory}G"

        script_content = textwrap.dedent(
            f"""\
            #!/usr/bin/env bash
            ##SBATCH --nodes=1                     # CURRENTLY UNSUPPORTED
            #SBATCH --ntasks-per-node={job.total_threads}       # number of tasks per node. Equivalent to ppn in TORQUE
            ##SBATCH --time=5-00:00:00                # CURRENTLY UNSUPPORTED
            #SBATCH --partition={job.queue}
            {job_memory_string}
            #SBATCH --job-name={self.settings.general.calculator}         # name of your job.
            
            set -e

            {scratch_partition_management_pre}{job.launch_command} > output.out 2> error.err
            {scratch_partition_management_post}
            """
        )

        return script_content


class PbsScheduler(SchedulerABC):
    def __init__(self, settings):
        super().__init__(settings)

    def execute(self):
        job_ids = []

        if not all([job.queue for job in self._pending_jobs]):
            # TODO: provide a more specific error pointing to which job is lacking a queue
            raise ValueError(
                "Jobs are missing a work queue."
                "Please make sure that all generators and calculators have a set work queue when using the PBS scheduler"
            )

        for job in self._pending_jobs:
            job.run(dry_run=True)
            job_dir = job.job_dir

            with open(os.path.join(job_dir, "submit_job.pbs"), "w") as pbs_file_handle:
                pbs_file_handle.write(self._pbs_script_content(job))

            job_id = subprocess.run(
                "qsub -terse submit_job.pbs'", cwd=job_dir, shell=True, capture_output=True, text=True
            )
            print(job_id)
            job_ids.append(job_id.stdout.strip("\n"))

        job_string = ":".join([f"{job_id}" for job_id in job_ids])

        print(job_string)
        dummy_job = subprocess.run(
            f"echo 'done' | qsub -W block=True depend=afterany:{job_string} --wrap='sleep 1'", shell=True
        )

        # Once completed, move the jobs out of the pending job list
        while self._pending_jobs:
            self._completed_jobs.append(self._pending_jobs.pop())

    def _pbs_script_content(self, job):

        scratch_dir = os.path.join(self.scratch_dir, "$USER")
        scratch_partition_management_pre = f"""
            # Create user-specific scratch directory if it does not exist
            if ! [ -d {scratch_dir} ]; then
            mkdir {scratch_dir}
            fi

            # Create a job-specific scratch directory and set environment variable
            export SCRDIR={scratch_dir}/$PBS_JOBID
            mkdir $SCRDIR

            ## Go to scratch directory and copy everything from working directory
            cd $SCRDIR
            cp -pr $PBS_O_WORKDIR/* . \n
            """

        scratch_partition_management_post = f"""
            # Copy results back to the original working directory
            cp -pr * $PBS_O_WORKDIR/

            # Cleanup scratch directory after successful job completion
            cd $PBS_O_WORKDIR
            rm -rf $SCRDIR
            """

        if self.scratch_dir is None:
            scratch_partition_management_pre = ""
            scratch_partition_management_post = ""

        if job.total_memory is None:
            job_memory_string = ""
        else:
            job_memory_string = f"#PBS -l mem={job.total_memory}gb"

        script_content = textwrap.dedent(
            f"""\
            #!/usr/bin/env bash
            #PBS -l nodes=1:ppn={job.total_threads}
            #SBATCH --ntasks-per-node={job.total_threads}       # number of tasks per node. Equivalent to ppn in TORQUE
            #PBS -q {job.queue}
            {job_memory_string}
            #PBS -N {self.settings.general.calculator}         # name of your job.
            
            set -e

            {scratch_partition_management_pre}{job.launch_command} > output.out 2> error.err
            {scratch_partition_management_post}
            """
        )

        return script_content
