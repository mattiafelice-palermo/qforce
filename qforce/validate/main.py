# Standard library imports
import os
from types import SimpleNamespace
import shutil

# Related third-party imports
import pkg_resources
from colt import Colt, from_commandline
from colt.validator import Validator

from .generators import AnnealerABC, _implemented_generators, _implemented_annealers, get_generator
from .calculators import Orca, get_calculator
from .archive import split_pdb_to_xyz
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

    # 2. Create scheduler
    scheduler = get_scheduler(parsed_settings)

    # 3. Create generator and run it through the scheduler
    generator = get_generator(parsed_settings)
    scheduler.add(generator)
    scheduler.execute()

    # 4. Create pool folder and move there the structures whose energy will be recalculated
    # TODO: at a certain point, this logic will be moved to a database manager
    if generator.structures_path is None:
        raise RuntimeError("No structures to resample has been found.")
    parsed_settings.general.pool_path = os.path.join(parsed_settings.general.job_dir, "pool")
    split_pdb_to_xyz(generator.structures_path, parsed_settings.general.pool_path)

    # 5. Create calculator and run resampling with it
    calculator = get_calculator(parsed_settings)
    scheduler.add(calculator)
    scheduler.execute()

    # 6. Correlate data and create plot

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
calculator = :: str

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

        # pprint(updated_settings)

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

        # CALCULATORS
        # print(questions["calculator"])
        # print({key: SimpleNamespace(**val) for key, val in questions.items()})

        questions.generate_block("orca", Orca.colt_user_input)

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

    return settings_file


def initialize_settings(config_file, presets=None):
    print("Remember to call the LOGO here")

    job_dir = os.getcwd()
    settings_file = _check_and_copy_settings_file(job_dir, config_file)

    settings = GeneralSettings.from_questions(config=settings_file, presets=presets, check_only=True)
    settings.general.job_dir = job_dir

    return settings
