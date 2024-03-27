from colt import Colt
from abc import ABC, abstractmethod
import shutil
import os
import subprocess
import textwrap
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
            # return XtbAnnealing(settings)
            raise NotImplementedError(f"Annealer '{annealer}' is not implemented.")
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

    #
    queue = :: str, optional
    """

    def __init__(self, settings):
        self.settings = settings
        self._t_min = settings.annealing.t_min
        self._t_max = settings.annealing.t_max
        self._n_steps = settings.annealing.n_steps
        self._dt = settings.annealing.dt
        self._n_configs = settings.annealing.n_configs
        self._job_type = "generator"
        self.queue = settings.annealing.queue
        self.generator_folder = self._create_working_folder()
        settings.general.number_of_structures = settings.annealing.n_configs

    def _create_working_folder(self):
        # Create generator working folder
        generator_folder = f"{self.settings.general.job_dir}/{self.settings.general.generator_method}"
        if os.path.exists(generator_folder):
            shutil.rmtree(generator_folder)
        os.makedirs(generator_folder)

        return generator_folder

    @abstractmethod
    def run(self): ...

    @property
    @abstractmethod
    def structures_path(self): ...

    def postprocess(self):
        pass

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

    # -1 uses all the available threads
    threads = -1 :: int, optional

    #
    total_memory = :: int, optional


    #
    queue = :: str, optional
    
    #
    custom_mdrun_flags = :: str, optional

    #
    custom_grompp_flags = :: str, optional

    #
    conda_environment = :: str, optional

    """

    def __init__(self, settings):
        super().__init__(settings)
        self.job_dir = self.generator_folder
        self.total_threads = self.settings.annealing.annealer.threads
        self._structures_path = None
        self.total_memory = settings.annealing.annealer.total_memory
        self.queue = settings.annealing.annealer.queue
        self.launch_command = "bash launch.sh"
        self.conda_environment = settings.annealing.annealer.conda_environment
        self.name = "GromacsAnnealing"
        # settings.general.number_of_structures += 1

        # Manage number of threads and custom gromacs flags
        if self.settings.annealing.annealer.threads == -1:
            threads = os.cpu_count()
            self.settings.annealing.annealer.threads = threads
            self.total_threads = threads
        else:
            threads = self.settings.annealing.annealer.threads

        self._setup_working_folder()

        self.mdp_annealing = textwrap.dedent(
            """; SIMULATED ANNEALING
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
        )

    # TODO: if not provided by user, use gmx insert-molecule to generate a box whose size is calculated from molecular size

    def _setup_working_folder(self):
        # Copy necessary files
        shutil.copy2(f"{self.settings.general.structure_file}", self.generator_folder)
        shutil.copy2(f"{self.settings.general.job_dir}/{self.settings.general.topology_file}", self.generator_folder)
        self._copy_mdp_file_to_generator_folder(self.generator_folder)

    def _copy_mdp_file_to_generator_folder(self, generator_folder: str) -> None:
        """Copy the MDP file to the generator folder."""
        user_mdp_filename = os.path.basename(self.settings.annealing.annealer.user_mdp_file)
        destination_path = os.path.join(generator_folder, "annealing.mdp")
        shutil.copy2(self.settings.annealing.annealer.user_mdp_file, destination_path)
        self.settings.annealing.annealer.system_mdp_file = destination_path

    def _generate_scripts(self) -> None:
        """Generate scripts to launch the calculation"""
        # Assign settings to local variables for increased readability
        script_path = os.path.join(self.generator_folder, "launch.sh")
        structure_file = self.settings.general.structure_file
        topology_file = self.settings.general.topology_file
        gromacs_executable = self.settings.annealing.annealer.gromacs_executable

        # Managing the pool folder creation here.. awkward, but it works
        pool_path = os.path.join(self.settings.general.job_dir, "pool")
        if os.path.exists(pool_path):
            shutil.rmtree(pool_path)
        os.makedirs(pool_path)

        if self.settings.annealing.annealer.custom_grompp_flags == None:
            custom_grompp_flags = ""
        else:
            custom_grompp_flags = self.settings.annealing.annealer.custom_grompp_flags

        if self.settings.annealing.annealer.custom_mdrun_flags == None:
            custom_mdrun_flags = ""
        else:
            custom_mdrun_flags = self.settings.annealing.annealer.custom_mdrun_flags

        box_vectors = self.settings.general.box_vectors.box_vectors

        # Write the composed script to the launch.sh file
        with open(script_path, "w") as script_handle:
            script_handle.write("#!/bin/sh\n\n")
            for directive in self.settings.annealing.annealer.custom_directives.split(";"):
                script_handle.write(directive + "\n")

            script_handle.write(
                f"gmx editconf -f {structure_file} -o molecule.gro -box {box_vectors[0]} {box_vectors[1]} {box_vectors[2]} > editconf.out 2> editconf.err\n"
            )

            script_handle.write(
                f"gmx grompp -f annealing.mdp -c molecule.gro -p {topology_file} -o annealing {custom_grompp_flags} > grompp.out 2> grompp.err\n"
            )

            script_handle.write(
                f"{gromacs_executable} mdrun -nt {self.total_threads} -deffnm annealing -c annealing {custom_mdrun_flags} > mdrun.out 2> mdrun.err\n"
            )

            # script_handle.write(f"echo 9 |  {gromacs_executable} energy -f annealing.edr -o annealing")

            script_handle.write(
                f"echo 0 | {gromacs_executable} trjconv -pbc whole -f annealing.trr -s annealing.tpr -o annealing.xtc > trjconv.out 2> trjconv.err\n"
            )

            script_handle.write(
                f"echo 0 0 | {gromacs_executable} trjconv -pbc mol -center -f annealing.xtc -s annealing.tpr -o annealing.pdb >> trjconv.out 2>> trjconv.err\n"
            )

    def _generate_input_files(self) -> None:
        """Generate input files by processing MDP settings and copying them to the specified folder.

        Args:
            generator_folder (str): The folder where the generated files will be stored.
        """
        self._replace_mdp_placeholders()  # fill annealing mdp section with user selected values
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

    def _update_system_mdp_file(self) -> None:
        """Update the system MDP file with correct settings."""
        trj_frequency = int(self.n_steps / self.n_configs)

        # -1 takes into account that the intial config is added to the trajectory (I think...)
        steps = trj_frequency * (self.n_configs - 1)

        # If user mistakenly specified gromacs keys related to annealing, remove them
        replacements = {
            "nsteps =": "nsteps".ljust(24) + "= " + f"{steps}",
            "dt =": "dt".ljust(24) + "= " + f"{self.dt}",
            "nstxout =": "nstxout".ljust(24) + "= " + f"{trj_frequency}",
            "nstenergy =": "nstenergy".ljust(24) + "= " + f"{trj_frequency}",
            "nstcalcenergy =": "1",
            "annealing =": "",
            "annealing_npoints =": "",
            "annealing_temp =": "",
            "annealing_time =": "",
            "gen-vel =": "",
            "gen-temp =": "",
            # "comm-mode =": "comm-mode = angular",
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

    def run(self, dry_run=False) -> None:
        self._generate_input_files()
        self._generate_scripts()

        if not dry_run:
            try:
                # TODO: launch command should be taken from the generator attribute
                mdrun = subprocess.Popen(
                    ["bash", "launch.sh"],
                    cwd=self.generator_folder,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )

                mdrun.wait()
            except Exception as e:
                print("Failed to run annealing run.\n", e)

        self.structures_path = os.path.join(self.generator_folder, "annealing.pdb")

    @property
    def structures_path(self):
        return self._structures_path

    @structures_path.setter
    def structures_path(self, argument):
        self._structures_path = argument


class XtbAnnealing(AnnealerABC):
    _user_input = """
    #
    test_input = sometest :: str
    """

    def generate_input_files(self):
        pass
