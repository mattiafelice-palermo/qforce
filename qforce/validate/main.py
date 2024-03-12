# Standard library imports
import os
from types import SimpleNamespace
import shutil

# Related third-party imports
import pkg_resources
from colt import Colt, from_commandline
from colt.validator import Validator

from .generators import AnnealerABC, _implemented_generators, _implemented_annealers, get_generator
from .calculators import get_calculators, get_calculator_class
from .archive import split_pdb_to_xyz
from .schedulers import SchedulerABC, get_scheduler
from .misc import get_fullpath, GroBoxEditor
from pprint import pprint

import re


def check_if_file_exists(filename):
    if not os.path.exists(filename) and not os.path.exists(f"{filename}_qforce"):
        raise ValueError(f'"{filename}" does not exist.\n')
    return filename


# define new Colt validator
Validator.overwrite_validator("file", check_if_file_exists)


@from_commandline(
    """
# File name for the optional settings.
config_file = :: file, optional, alias=o
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
def run_validator(config_file):
    # 1. Read settings from the user provided file
    settings = initialize_settings(config_file)

    # 2. Create scheduler
    scheduler = get_scheduler(settings)

    # 3. Create generator and run it through the scheduler
    # TODO: allow for multiple generators run
    generator = get_generator(settings)
    scheduler.add(generator)
    scheduler.execute()

    # 4. Create pool folder and move there the structures whose energy will be recalculated
    # TODO: at a certain point, this logic will be moved to a database manager
    if generator.structures_path is None:
        raise RuntimeError("No structures to resample has been found.")
    settings.general.pool_path = os.path.join(settings.general.job_dir, "pool")
    split_pdb_to_xyz(generator.structures_path, settings.general.pool_path)

    # 5. Create calculators and run resamplings
    calculators = get_calculators(settings)
    for calculator in calculators:
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
calculators = :: str

#
scheduler = :: str :: [manual, system, pbs, slurm]

#
shell = bash :: str, optional :: [bash, sh, zsh]
"""

    @classmethod
    def from_config(cls, settings):
        settings_from_questions = cls._set_config(settings)
        reformatted_settings = cls.reformat_settings(settings_from_questions)
        return reformatted_settings

    @staticmethod
    def _set_config(settings):
        # NOTE: Data from colt.answers can be explored with settings[value].to_dict()
        # Hacking colt to restructure calculators section - colt does not allow to directly edit answers
        # so I access the underlying class dictionary to access the "_data" attribute, which is a dictionary itself
        calculator_jobs_string = settings["general"]["calculators"]
        calculator_jobs = parse_jobs_to_dict(calculator_jobs_string)

        settings.__dict__["_data"]["calculators"] = {}

        for job_id, calculator_name in calculator_jobs.items():
            job_string = f"{job_id}::{calculator_name}"
            print(job_string)
            calculator_settings = SimpleNamespace(**settings[job_string])
            settings.__dict__["_data"]["calculators"].update({job_string: calculator_settings})
            del settings.__dict__["_data"][job_string]

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

    def reformat_settings(settings_from_questions):
        pprint(settings_from_questions.calculators)

        for calculator_name, calculator_settings in vars(settings_from_questions.calculators).items():
            calculator_settings.name = calculator_name.split("::")[0]
            calculator_settings.driver = calculator_name.split("::")[1]

        # Build absolute path for the structure file
        settings_from_questions.general.structure_file = get_fullpath(settings_from_questions.general.structure_file)

        settings_from_questions.general.box_vectors = GroBoxEditor(settings_from_questions.general.structure_file)
        return settings_from_questions

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

        calculator_jobs = {}

        with open(cls._config_file, "r") as config_handle:
            for line in config_handle:
                if "calculator" in line:
                    calculator_jobs = parse_jobs_to_dict(line)

        print(calculator_jobs)

        # CALCULATORS
        for job_id, calculator_name in calculator_jobs.items():
            calc = get_calculator_class(calculator_name)
            questions.generate_block(f"{job_id}::{calculator_name}", calc.colt_user_input)

        # questions.generate_block("scan", DihedralScan.colt_user_input)
        # questions.generate_cases(
        #     "annealing",
        #     {key: software.colt_user_input for key, software in implemented_qm_software.items()},
        #     block="qm",
        # )
        # questions.generate_block("terms", Terms.get_questions())

    @classmethod
    def from_questions(
        cls, *args, check_only=False, ask_all=False, ask_defaults=True, config=None, presets=None, **kwargs
    ):
        cls._config_file = config
        out = super().from_questions(
            *args,
            check_only=check_only,
            ask_all=ask_all,
            ask_defaults=ask_defaults,
            config=config,
            presets=presets,
            **kwargs,
        )
        return out


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


def parse_jobs_to_dict(line):
    # Verify the overall format is correct (jobs separated by commas, spaces allowed)
    clean_line = re.sub(r"^\s*calculators\s*=\s*", "", line, flags=re.IGNORECASE)

    if not re.match(r"^\s*(\w+::\w+\s*,\s*)*\w+::\w+\s*$", clean_line):
        raise ValueError("The line is not formatted correctly. Expected format: 'key::value, key2::value2, ...'")

    # Pattern to extract key-value pairs, ignoring spaces around ::
    pattern = r"\s*(\w+)::(\w+)\s*"
    matches = re.findall(pattern, line)

    if matches:
        return dict(matches)
    else:
        raise ValueError("No valid 'key::value' pairs found.")
