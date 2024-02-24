from colt import Colt
from abc import ABC, abstractmethod
import shutil
import os

from pprint import pprint

_implemented_generators = ["annealing", "crest", "qcg"]


# Used in main.py to read the available annealers
def _implemented_annealers():
    annealers = {"gromacs": GromacsAnnealing, "xtb": XtbAnnealing}
    return annealers.items()


def get_generator(settings):
    """
    Returns the appropriate generator based on the specified settings.

    Args:
        settings: The settings object containing the configurations for qforce-validate

    Returns:
        Generator: An instance of a generator based on the specified method.

    Raises:
        NotImplementedError: If the generator method specified in the settings is not supported.
    """
    generator_method = settings.general.generator_method.lower()

    if generator_method == "annealing":
        annealer = settings.annealing.annealer.name.lower()

        if annealer == "gromacs":
            return GromacsAnnealing(settings)
        elif annealer == "xtb":
            return XtbAnnealing(settings)
        else:
            raise NotImplementedError(f"Annealer '{annealer}' is not implemented.")

    raise NotImplementedError(f"Generator method '{generator_method}' is not implemented.")


class AnnealerABC(ABC, Colt):
    _user_input = """
    # 
    t_min = 200. :: float

    # 
    t_max = 500. :: float
    
    # number of steps
    n_steps = 10000000 :: int

    # dt (ps)
    dt = 0.001 :: float

    # number of configurations
    n_configs = 1000 :: int
    """

    def __init__(self, settings):
        self.settings = settings
        self._t_min = settings.annealing.t_min
        self._t_max = settings.annealing.t_max
        self._n_steps = settings.annealing.n_steps
        self._dt = settings.annealing.dt
        self._n_configs = settings.annealing.n_configs
        self._job_type = "generator"

    @abstractmethod
    def generate_input_files(self): ...

    @property
    def t_min(self):
        return self._t_min

    @property
    def t_max(self):
        return self._t_max

    @property
    def n_temps(self):
        return self._n_temps

    @property
    def n_steps(self):
        return self._n_steps

    @property
    def dt(self):
        return self._dt

    @property
    def n_configs(self):
        return self._n_configs

    @property
    def job_type(self):
        return self._job_type


class GromacsAnnealing(AnnealerABC):
    _user_input = """
    # 
    user_mdp_file = annealing.mdp :: str, optional
    gromacs_executable = gmx :: str, optional
    
    # e.g. module load gromacs; source /usr/local/gromacs/bin/GMXRC or others
    custom_directives = module load gromacs :: str, optional
    """

    def __init__(self, settings):
        super().__init__(settings)
        self.generator_folder = f"{settings.general.job_dir}/{settings.general.generator_method}"
        self.mdp_annealing = """; SIMULATED ANNEALING
; Type of annealing for each temperature group (no/single/periodic)
annealing                = single
; Number of time points to use for specifying annealing in each group
annealing_npoints        = 2
; List of times at the annealing points for each group
annealing_time           = 0 ANNEALING_TIME
; Temp. at each annealing point, for each group.
annealing_temp           = ANNEALING_TEMP

; GENERATE VELOCITIES FOR STARTUP RUN
gen-vel                  = yes
gen-temp                 = INITIAL_TEMP"""

    # TODO: if not provided by user, use gmx insert-molecule to generate a box whose size is calculated from molecular size

    def generate_scripts(self) -> None:
        """Generate scripts to launch the calculation"""
        script_path = os.path.join(self.generator_folder, "launch.sh")
        job_dir = self.settings.general.job_dir
        structure_file = self.settings.general.structure_file
        topology_file = self.settings.general.topology_file
        system_mdp_file = self.settings.annealing.annealer.system_mdp_file
        gromacs_executable = self.settings.annealing.annealer.gromacs_executable
        with open(script_path, "w") as script_handle:
            script_handle.write("#!/bin/sh\n\n")
            for directive in self.settings.annealing.annealer.custom_directives.split(";"):
                script_handle.write(directive + "\n")
            script_handle.write(
                f"gmx grompp -f {system_mdp_file} -c {job_dir}/{structure_file} -p {job_dir}/{topology_file} -o annealing &> grompp.out\n"
            )

            script_handle.write(
                f"{gromacs_executable} mdrun -nt 4 -deffnm annealing -c annealing -tunepme  &> mdrun.out &\n"
            )

    def generate_input_files(self) -> None:
        """Generate input files by processing MDP settings and copying them to the specified folder.

        Args:
            generator_folder (str): The folder where the generated files will be stored.
        """
        self._replace_mdp_placeholders()  # fill annealing mdp section with user selected values
        self._copy_mdp_file_to_generator_folder(self.generator_folder)
        self._update_system_mdp_file()  # create the final mdp file for annealing, superseeding some user input

    def _replace_mdp_placeholders(self) -> None:
        """Replace placeholders in MDP string with actual values."""
        # Calculate annealing time
        annealing_time = self.n_steps * self.dt

        # Replace placeholders in MDP string with actual values
        replacements = {
            "ANNEALING_TIME": str(annealing_time),
            "ANNEALING_TEMP": f"{self.t_min} {self.t_max}",
            "INITIAL_TEMP": str(self.t_min),
        }

        for key, value in replacements.items():
            self.mdp_annealing = self.mdp_annealing.replace(key, value)

    def _copy_mdp_file_to_generator_folder(self, generator_folder: str) -> None:
        """Copy the MDP file to the generator folder."""
        user_mdp_filename = os.path.basename(self.settings.annealing.annealer.user_mdp_file)
        print(user_mdp_filename, generator_folder)
        destination_path = os.path.join(generator_folder, user_mdp_filename)
        shutil.copy2(self.settings.annealing.annealer.user_mdp_file, destination_path)
        self.settings.annealing.annealer.system_mdp_file = destination_path

    def _update_system_mdp_file(self) -> None:
        """Update the system MDP file with correct settings."""
        trj_frequency = int(self.n_steps / self.n_configs)
        # If user mistakenly specified gromacs keys related to annealing, remove them
        replacements = {
            "nsteps =": "nsteps".ljust(24) + "= " + f"{self.n_steps}",
            "dt =": "dt".ljust(24) + "= " + f"{self.dt}",
            "nstxout =": "nstxout".ljust(24) + "= " + f"{trj_frequency}",
            "annealing =": "",
            "annealing_npoints =": "",
            "annealing_temp =": "",
            "annealing_time =": "",
            "gen-vel =": "",
            "gen-temp =": "",
        }

        system_mdp_lines = []

        # find instances of mdp keys to replace in user input and replace them
        with open(self.settings.annealing.annealer.system_mdp_file, "r") as system_mdp_file:
            found = []  # keep track of the mdp keys that have been found in the user input
            for line in system_mdp_file:
                line = " ".join(line.split())  # remove spurious whitespaces
                for key, value in replacements.items():
                    # only consider content before ;, in weird cases a key is found in a comment...
                    splitted = line.split(";", 1)[0]
                    if (key in splitted) or (key.replace(" ", "") in splitted):
                        system_mdp_lines.append(value + "\n")
                        found.append(key)
                        break  # keys must appear only once in the mdp file, no need to keep checking
                else:
                    system_mdp_lines.append(line + "\n")

        # If a keys to be replaced are not found in the user input, add them
        not_found = set(replacements) - set(found)

        for key in not_found:
            system_mdp_lines.append(replacements[key] + "\n")

        # add annealing section to mdp file
        for line in self.mdp_annealing.split("\n"):
            system_mdp_lines.append(line + "\n")

        with open(self.settings.annealing.annealer.system_mdp_file, "w") as file:
            file.writelines(system_mdp_lines)


class XtbAnnealing(AnnealerABC):
    _user_input = """
    #
    test_input = sometest :: str
    """

    def generate_input_files(self):
        pass
