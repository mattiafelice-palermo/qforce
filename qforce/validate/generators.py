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


class GromacsAnnealing(AnnealerABC):
    _user_input = """
    # 
    user_mdp_file = annealing.mdp :: str, optional
    """

    def __init__(self, settings):
        super().__init__(settings)
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

    def generate_input_files(self, generator_folder: str) -> None:
        """Generate input files by processing MDP settings and copying them to the specified folder.

        Args:
            generator_folder (str): The folder where the generated files will be stored.
        """
        # Calculate annealing time
        annealing_time = self.n_steps * self.dt
        trj_frequency = self.n_steps / self.n_configs

        # Replace placeholders in MDP string with actual values
        replacements = {
            "ANNEALING_TIME": str(annealing_time),
            "ANNEALING_TEMP": f"{self.t_min} {self.t_max}",
            "INITIAL_TEMP": str(self.t_min),
        }

        for key, value in replacements.items():
            self.mdp_annealing = self.mdp_annealing.replace(key, value)

        # Copy the MDP file to the generator folder
        user_mdp_filename = os.path.basename(self.settings.annealing.annealer.user_mdp_file)
        destination_path = os.path.join(generator_folder, user_mdp_filename)
        shutil.copy2(self.settings.annealing.annealer.user_mdp_file, destination_path)
        self.settings.annealing.annealer.system_mdp_file = destination_path

        # If user mistakenly specified gromacs keys related to annealing, remove them
        replacements = {
            "nsteps =": "nsteps".ljust(24) + "= " + f"{self.n_steps}",
            "dt =": "dt".ljust(24) + "= " + f"{self.dt}",
            "nstxout =": "nstxout".ljust(24) + "= " + f"{trj_frequency}",
        }

        system_mdp_lines = []
        with open(self.settings.annealing.annealer.system_mdp_file, "r") as system_mdp_file:
            found = []
            for line in system_mdp_file:
                line = " ".join(line.split())  # remove spurious whitespaces
                for key, value in replacements.items():
                    splitted = line.split(";", 1)

                    if (key in splitted[0]) or (key.replace(" ", "") in splitted[0]):
                        # if the line that matches the key starts with ;, skip it!
                        system_mdp_lines.append(value)
                        found.append(key)
                        break  # keys must appear only once in the mdp file, no need to keep checking
                else:
                    system_mdp_lines.append(line)

        not_found = set(replacements.keys()).difference(found)

        for key in not_found:
            system_mdp_lines.append(replacements[key])

        for line in self.mdp_annealing.split("\n"):
            system_mdp_lines.append(line)


class XtbAnnealing(AnnealerABC):
    _user_input = """
    #
    test_input = sometest :: str
    """

    def generate_input_files(self):
        print(self._t_min)
