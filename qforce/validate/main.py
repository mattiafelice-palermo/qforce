# Standard library imports
import os
from types import SimpleNamespace

# Related third-party imports
import pkg_resources
from colt import Colt, from_commandline
from colt.validator import Validator


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
    # read settings from the user provided file
    parsed_settings = initialize_settings(settings)
    # print(parsed_settings)
    # 1. Read settings

    # 2. Launch generator
    #  - if generator is based on gromacs, energies should be stored

    # 3. Energy sorting and duplicates removal

    # 4. Resampling

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
generator_method = :: str :: [anneal, anneal_xtb, crest, qcg_microsolv]

#
scheduler = :: str :: [none, manual, pbs, slurm]
"""

    @staticmethod
    def _set_config(settings):
        settings.update({key: SimpleNamespace(**val) for key, val in settings.items()})
        return SimpleNamespace(**settings)

    @classmethod
    def from_config(cls, settings):
        return cls._set_config(settings)

    # @classmethod
    # def _extend_user_input(cls, questions):
    #     questions.generate_block("qm", QM.colt_user_input)
    #     questions.generate_block("scan", DihedralScan.colt_user_input)
    #     questions.generate_cases(
    #         "software",
    #         {key: software.colt_user_input for key, software in implemented_qm_software.items()},
    #         block="qm",
    #     )
    #     questions.generate_block("terms", Terms.get_questions())


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

    settings_file = os.path.join(job_dir, "settings.ini")

    # if config_file is not None:
    #     if isinstance(config_file, StringIO):
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

    settings = GeneralSettings.from_questions(
        config=settings_file, presets=presets, check_only=True
    )

    return settings
