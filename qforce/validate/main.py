# Standard library imports
import os
from types import SimpleNamespace
import shutil

# Related third-party imports
import pkg_resources
from colt import Colt, from_commandline
from colt.validator import Validator

from .generators import AnnealerABC, _implemented_generators, _implemented_annealers, get_generator
from .schedulers import SchedulerABC, get_scheduler
from pprint import pprint


def check_if_file_exists(filename):
    if not os.path.exists(filename) and not os.path.exists(f"{filename}_qforce"):
        raise ValueError(f'"{filename}" does not exist.\n')
    return filename


# define new Colt validator
Validator.overwrite_validator("file", check_if_file_exists)


@from_commandline(
    """
# File name for the optional settings.
settings = :: file, optional, alias=o
""",
    description={
        "logo": "Add qforce LOGO here",
        "alias": "qforce",
        "arg_format": {
            "name": 12,
            "comment": 60,
        },
    },
)

# Entry point for the qforce-validate routine
def run_validator(settings):
    # 1. Read settings from the user provided file
    parsed_settings = initialize_settings(settings)

    # Create generator working folder
    generator_folder = f"{parsed_settings.general.job_dir}/{parsed_settings.general.generator_method}"
    if os.path.exists(generator_folder):
        shutil.rmtree(generator_folder)
    os.makedirs(generator_folder)

    generator = get_generator(parsed_settings)
    scheduler = get_scheduler(parsed_settings)
    scheduler.add(generator)
    scheduler.execute()

    # Generate generator input and launching scripts
    if scheduler == "manual":
        pass
    if scheduler == "auto":  # maybe system is better?
        # 1. generate generator inputs and scripts
        # 2. launch locally using available resources
        pass
    if scheduler == "pbs":
        # 1. generate generator inputs and scripts
        # 2. generate generator inputs and scripts
        # 3. generate queue files
        # 4. submit to queue
        # 5. store queue id for dependeny
        pass

    # 2.2 scheduler = none
    #  - System settings for MD and QM calculation should be provided (tasks, memory per task)

    # 2.3 scheduler = pbs
    #  - Scheduler settings for MD and QM calculation should be provided (tasks, memory per task, possibly nodes, queue names)

    # 2.4 scheduler = slurm
    #  - Scheduler settings for MD and QM calculation should be provided (tasks, memory per task, possibly nodes, queue names)

    # 3. Run generator
    #  - Run differently based on the scheduler option

    # 3. Energy sorting and duplicates removal

    # 4. Resampling
    #  - Two resampling techniques: GROMACS MD and any QM software
    #  - Simulation settings for MD and QM calculation should be provided
    #  - Scheduler settings for MD and QM calculation should be provided (tasks, memory per task, possibly nodes, queue names)

    pass


class GeneralSettings(Colt):
    _user_input = """
[general]
# 
source = auto :: str :: [auto, user]

#
structure_file = :: str

#
topology_file = :: str

#
md_settings_file = :: str

#
generator_method = :: str :: [annealing, crest, qcg_microsolv]

#
scheduler = :: str :: [manual, system, pbs, slurm]
"""

    @staticmethod
    def _set_config(settings):
        # Extracting annealer settings and updating its format
        annealer_info = settings["annealing"]["annealer"]
        annealer_name = annealer_info.value
        annealer_dict = dict(annealer_info)  # Convert Colt object to dict
        annealer_dict["name"] = annealer_name

        # Converting annealer settings to SimpleNamespace for easier access
        annealer_settings = SimpleNamespace(**annealer_dict)
        settings["annealing"].update({"annealer": annealer_settings})  # update colt.questions object

        # Converting all settings to SimpleNamespace for uniformity and ease of use
        updated_settings = {key: SimpleNamespace(**val) for key, val in settings.items()}

        return SimpleNamespace(**updated_settings)

    @classmethod
    def from_config(cls, settings):
        return cls._set_config(settings)

    @classmethod
    def _extend_user_input(cls, questions):
        # GENERATORS
        questions.generate_block("annealing", AnnealerABC.colt_user_input)
        annealers = {
            annealer_name: annealer_class.colt_user_input for annealer_name, annealer_class in _implemented_annealers()
        }  # colt_user_input returns the parameters set in the annealer class _user_input class variable e.g. GromacsAnnealer

        questions.generate_cases("annealer", annealers, block="annealing")

        # SCHEDULERS
        questions.generate_block("scheduler", SchedulerABC.colt_user_input)

        # questions.generate_block("scan", DihedralScan.colt_user_input)
        # questions.generate_cases(
        #     "annealing",
        #     {key: software.colt_user_input for key, software in implemented_qm_software.items()},
        #     block="qm",
        # )
        # questions.generate_block("terms", Terms.get_questions())


def _get_job_info(filename):
    job = {}
    filename = filename.rstrip("/")
    base = os.path.basename(filename)
    path = os.path.dirname(filename)
    if path != "":
        path = f"{path}/"

    if os.path.isfile(filename):
        job["coord_file"] = filename
        job["name"] = base.split(".")[0]
    else:
        job["coord_file"] = False
        job["name"] = base.split("_qforce")[0]

    job["dir"] = f'{path}{job["name"]}_qforce'
    job["frag_dir"] = f'{job["dir"]}/fragments'
    job["md_data"] = pkg_resources.resource_filename("qforce", "data")
    os.makedirs(job["dir"], exist_ok=True)
    return SimpleNamespace(**job)


def _check_and_copy_settings_file(job_dir, config_file):
    """
    If options are provided as a file, copy that to job directory.
    If options are provided as StringIO, write that to job directory.
    """

    settings_file = os.path.join(job_dir, ".settings.parsed.ini")
    shutil.copy2(config_file, settings_file)

    # if config_file is not None:
    #     if isinstance(config_file, StringIO):q
    #         with open(settings_file, "w") as fh:
    #             config_file.seek(0)
    #             fh.write(config_file.read())
    #     else:
    #         shutil.copy2(config_file, settings_file)

    return settings_file


def initialize_settings(config_file, presets=None):
    print("Remember to call the LOGO here")

    job_dir = os.getcwd()
    settings_file = _check_and_copy_settings_file(job_dir, config_file)

    settings = GeneralSettings.from_questions(config=settings_file, presets=presets, check_only=True)
    settings.general.job_dir = job_dir

    return settings
