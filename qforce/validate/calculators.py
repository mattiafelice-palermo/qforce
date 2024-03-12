from colt import Colt
from .misc import get_fullpath
from .archive import extract_energies

from abc import ABC, abstractmethod
import shutil
import os
import textwrap
import subprocess
import inspect


def get_calculators(settings):
    """
    Returns the appropriate calculator based on the specified settings.

    Args:
        settings: The settings object containing the configurations for qforce-validate

    Returns:
        Calculator: An instance of a calculator based on the specified options.

    Raises:
        NotImplementedError: If the calculator specified in the settings is not supported.
    """
    calculators = []

    # Create calculator objects from the user input and return them in a list
    for job_string, calculator_settings in vars(settings.calculators).items():  # fetch multiple jobs
        if calculator_settings.driver == "orca":
            calculators.append(Orca(job_string, settings))
        elif calculator_settings.driver == "gromacs":
            calculators.append(Gromacs(job_string, settings))
        else:  # in theory we should never get here, but just to be safe...
            raise NotImplementedError(f"Calculator '{calculator_settings.driver}' is not implemented.")

    return calculators


def get_calculator_class(name):
    # Convert the user input to lowercase for case-insensitive comparison
    name_lowercase = name.lower()

    # Filter the globals for subclasses of CalculatorABC
    # and prepare a dictionary mapping lowercase class names to their actual class objects
    calculator_classes = {
        cls_name.lower(): cls_obj
        for cls_name, cls_obj in globals().items()
        if inspect.isclass(cls_obj) and issubclass(cls_obj, CalculatorABC) and cls_obj is not CalculatorABC
    }

    # Try to find a class that matches the user input (case-insensitive)
    for cls_name_lower, cls_obj in calculator_classes.items():
        if cls_name_lower == name_lowercase:
            return cls_obj

    # If no class was found that matches the user input, raise an error
    raise NotImplementedError(f"Calculator '{name}' is not implemented.")


class CalculatorABC(ABC):
    def __init__(self, settings):
        self.settings = settings
        self._check_files_in_pool()

    def _check_files_in_pool(self):
        pool_path = self.settings.general.pool_path
        for file in os.listdir(pool_path):
            if file.endswith(".xyz"):
                return

        raise FileNotFoundError("No xyz file has been found in the pool directory.")

    @abstractmethod
    def run(self):
        raise NotImplementedError

    @abstractmethod
    def postprocess(self):
        raise NotImplementedError


class Orca(CalculatorABC, Colt):
    """Just instruct generators and resamplers to write input files and launch scripts to disk."""

    _user_input = """
    #
    charge = 0 :: str

    #
    multiplicity = 1 :: str

    #
    method = B97-3c :: str

    #
    blocks = :: str, optional

    #
    total_threads = :: int

    # In gygabytes
    total_memory = :: int

    #
    single_calculation_threads = 1 :: int

    #
    orca_executable = orca :: str, optional

    # list of modules to be loaded (semicolon separated)
    modules = :: str, optional

    # list of path variable and path to be prepended in the format SOMEVARIABLEPATH: SOMEPATH (semicolon separated)
    # example: PATH: /usr/lib64/openmpi/bin; LD_LIBRARY_PATH: /usr/lib64/openmpi/lib
    exports = :: str, optional

    #
    queue = :: str, optional

    #
    conda_environment = ::str, optional

    """

    def __init__(self, job_string, settings):
        super().__init__(settings)

        # Automatically set all attributes from calculator_settings
        calculator_settings = vars(settings.calculators)[job_string]

        for attr, value in vars(calculator_settings).items():
            setattr(self, attr, value)

        # Other calculator settings that do not come directly from the colt calculator settings
        self.calculator_folder = os.path.join(self.settings.general.job_dir, job_string).replace("::", "__")
        self.launch_command = "python dispatcher.py"
        self.job_dir = self.calculator_folder  # just for the schedulers

    def run(self, dry_run=False):
        self._setup_working_folder()
        self._generate_orca_template()
        self._generate_launch_script()

        # For the system scheduler
        if not dry_run:
            try:
                run_orca = subprocess.Popen(
                    ["python", "dispatcher.py"],
                    cwd=self.calculator_folder,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )

                run_orca.wait()
            except Exception as e:
                print("Failed to run Orca calculator.\n", e)

    def _setup_working_folder(self):
        # Create generator working folder
        if os.path.exists(self.calculator_folder):
            shutil.rmtree(self.calculator_folder)
        os.makedirs(self.calculator_folder)

        # TODO: generate also the nth folders for the calculation? they will need to be copied though...

    def _generate_launch_script(self):
        calculator_folder = self.calculator_folder
        # Create the job script based on user input and template
        script_content = self._python_script_content()
        with open(os.path.join(calculator_folder, "dispatcher.py"), "w") as file:
            file.write(script_content)

    def _python_script_content(self):
        """
        Generates a Python script content for ORCA job execution and management.

        This method constructs a multi-line string that outlines a Python script. The script is designed to
        dynamically create directories for ORCA calculations, manage input and output files, execute ORCA
        with predefined settings, and aggregate results. The script leverages a job dispatcher for parallel
        execution of ORCA jobs based on the number of molecular structures defined in the settings.

        Returns:
            str: A dedented multi-line string containing the Python script for executing ORCA jobs,
                 handling file operations, and aggregating results. The script includes functions for
                 directory creation, ORCA execution setup, and output management, configured according
                 to instance attributes and settings.

        The generated script uses external tools (`sed`, `orca`) and assumes their availability in the execution environment.
        It dynamically adjusts to the number of structures and computational resources specified in the settings.
        """

        # This is where the script is customized based on the user input
        number_of_structures = self.settings.general.number_of_structures
        modules_string = build_modules_string(self.modules)
        exports_string = build_exports_string(self.exports)

        return textwrap.dedent(
            f"""
        import jobdispatcher as jd
        import subprocess
        import os
        import shutil

        def run_orca(filename):
            # Call the orca_calc.sh script and capture its output
            folder_path = create_directory(os.path.join("{self.calculator_folder}", filename))
            shutil.copy2('orca_template.inp', os.path.join(folder_path, f"{{filename}}.inp"))
            xyz_path = os.path.join("{self.settings.general.pool_path}", f"{{filename}}.xyz")
            subprocess.run(f"sed -i 's#FILENAME#{{xyz_path}}#g' {{filename}}.inp", shell=True, cwd=folder_path, check=True)

            result = subprocess.run(f'{modules_string}{exports_string} $(which orca) {{filename}}.inp > {{filename}}.out 2> {{filename}}.err', shell=True, cwd=folder_path, check=True)

            print(result)

            return result.stdout

        def create_directory(folder):
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder)
            return folder

        jobs = []
        for i in range(1, {number_of_structures}+1):
            job = jd.Job(name=f"molecule_{{i}}", function=lambda i=i: run_orca(f"molecule_{{i}}"), cores={self.single_calculation_threads})
            jobs.append(job)

        dispatcher = jd.JobDispatcher(jobs, maxcores={self.total_threads}, engine="multiprocessing")
        results = dispatcher.run()

        # Writing the outputs to a file
        with open("orca_outputs.txt", "w") as outfile:
            for result in results.values():
                outfile.write(result)
        """
        )

    def _generate_orca_template(self):
        calculator_folder = self.calculator_folder

        orca_input_path = os.path.join(calculator_folder, "orca_template.inp")

        if self.single_calculation_threads > 1:
            pal_string = f"PAL{self.single_calculation_threads}"
        else:
            pal_string = ""

        template_string = textwrap.dedent(
            f"""
        ! {self.method} {pal_string}
        * xyzfile {self.charge} {self.multiplicity} FILENAME
        """
        )

        with open(orca_input_path, "w") as orca_input_handle:
            orca_input_handle.write(template_string)

    def postprocess(self):
        energies = {}
        for molecule_index in range(1, self.settings.general.number_of_structures + 1):
            with open(
                os.path.join(self.calculator_folder, f"molecule_{molecule_index}/molecule_{molecule_index}.out"), "r"
            ) as file_handle:
                for line in file_handle:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        energies[molecule_index] = float(line.split()[4])

        print(energies)
        return energies


class Gromacs(CalculatorABC, Colt):
    """Just instruct generators and resamplers to write input files and launch scripts to disk."""

    _user_input = """
    #
    topology_file = :: str

    #
    md_settings_file = :: str

    #
    total_threads = :: int

    # In gygabytes
    total_memory = :: int

    #
    custom_mdrun_flags = :: str, optional

    #
    custom_grompp_flags = :: str, optional

    #
    conda_environment = :: str, optional
    
    # The number of structures per thread is calculated based on single_calculation_threads and total_threads
    single_calculation_threads = 1 :: int

    # 
    gromacs_executable = gmx :: str, optional

    # list of modules to be loaded (semicolon separated)
    modules = :: str, optional

    # list of path variable and path to be prepended in the format SOMEVARIABLEPATH: SOMEPATH (semicolon separated)
    # example: PATH: /usr/lib64/openmpi/bin; LD_LIBRARY_PATH: /usr/lib64/openmpi/lib
    exports = :: str, optional

    #
    queue = :: str, optional

    #
    conda_environment = ::str, optional

    """

    def __init__(self, job_string, settings):
        super().__init__(settings)

        # Automatically set all attributes from calculator_settings
        calculator_settings = vars(settings.calculators)[job_string]

        for attr, value in vars(calculator_settings).items():
            setattr(self, attr, value)

        self.topology_file = get_fullpath(self.topology_file)
        self.md_settings_file = get_fullpath(self.md_settings_file)

        # Other calculator settings that do not come directly from the colt calculator settings
        self.calculator_folder = os.path.join(self.settings.general.job_dir, job_string).replace("::", "__")
        self.launch_command = "python dispatcher.py"
        self.job_dir = self.calculator_folder  # just for the schedulers
        self.number_of_batches = self.total_threads // self.single_calculation_threads
        self.structures_per_batch = self.settings.general.number_of_structures // self.number_of_batches

    def run(self, dry_run=False):
        self._setup_working_folder()
        self._generate_launch_script()

        # For the system scheduler
        if not dry_run:
            try:
                command = "python dispatcher.py > dispatcher.out 2> dispatcher.err"
                subprocess.run(command, shell=True, cwd=self.calculator_folder, check=True, executable="/bin/bash")
            except Exception as e:
                print("Failed to run GROMACS calculator.\n", e)

    def _setup_working_folder(self):
        # Create generator working folder
        if os.path.exists(self.calculator_folder):
            shutil.rmtree(self.calculator_folder)
        os.makedirs(self.calculator_folder)

    def _generate_launch_script(self):
        # Create the job script based on user input and template
        script_content = self._python_script_content()
        with open(os.path.join(self.calculator_folder, "dispatcher.py"), "w") as file:
            file.write(script_content)

    def _python_script_content(self):
        # This is where the script is customized based on the user input
        number_of_batches = self.number_of_batches
        structures_per_batch = self.structures_per_batch
        modules_string = build_modules_string(self.modules)
        exports_string = build_exports_string(self.exports)
        box_vectors = self.settings.general.box_vectors.get_original_box_vector_string()

        return textwrap.dedent(
            f"""
        import jobdispatcher as jd
        from openbabel import pybel
        import subprocess
        import os
        import shutil

        def run_gromacs(batch_number):
            # Call the orca_calc.sh script and capture its output
            folder_path = create_directory(os.path.join("{self.calculator_folder}", f"batch_{{batch_number}}"))
            batchfile = generate_gro_file(batch_number, folder_path)
 
            # Copy necessary files into the running folder (maybe not necessary, just refer to root folder?)

            # Run GROMPP
            try: 
                grompp = subprocess.run(f'{modules_string}{exports_string} $(which {self.gromacs_executable}) grompp -f resample.mdp'
                                        ' -c resample.gro -p resample.top -o resample.tpr -maxwarn 2 > grompp.out 2> grompp.err', 
                                        shell=True, cwd=folder_path, check=True, executable='/bin/bash')
                print(grompp)
            except:
                raise RuntimeError(f"Failed to run GROMACS grompp command in {self.name}::{self.driver} for batch {{batch_number}}."
                                    "Check the calculator folder to inspect the error:"
                                    " -> {self.calculator_folder}")

            # Run MDRUN
            try:
                mdrun = subprocess.run(f'{modules_string}{exports_string} $(which {self.gromacs_executable}) mdrun -nt {self.single_calculation_threads}'
                                       f' -s resample.tpr -rerun {{batchfile}} > mdrun.out 2> mdrun.err', shell=True, cwd=folder_path, check=True, executable='/bin/bash')

                print(mdrun)
            except:
                raise RuntimeError(f"Failed to run GROMACS mdrun command in {self.name}::{self.driver} for batch {{batch_number}}."
                                    "Check the calculator folder to inspect the error:"
                                    " -> {self.calculator_folder}")
            return #result.stdout
        
        def generate_gro_file(batch_number, folder_path):
            initial_molecule_id = 1+{structures_per_batch}*(batch_number-1)
            final_molecule_id = initial_molecule_id + {structures_per_batch}
            batchfile_path = os.path.join(folder_path, f"batch_{{batch_number}}.gro")

            with open(batchfile_path, 'w') as gro_file:
                for i in range(initial_molecule_id, final_molecule_id):
                    mol = next(pybel.readfile('xyz', os.path.join("{self.settings.general.pool_path}", f"molecule_{{i}}.xyz")))
                    gro_str = mol.write('gro')
                    gro_lines = gro_str.splitlines()
                    gro_lines[-1] = "{box_vectors}"
                    modified_gro_str = "\\n".join(gro_lines)
                    gro_file.write(modified_gro_str + "\\n")  # Add a newline to separate entries
            return batchfile_path



        def create_directory(folder):
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder)
            os.symlink("{self.topology_file}", os.path.join(folder, "resample.top"))
            os.symlink("{self.md_settings_file}", os.path.join(folder, "resample.mdp"))
            os.symlink("{self.settings.general.structure_file}", os.path.join(folder, "resample.gro"))
            return folder

        jobs = []
        for i in range(1, {number_of_batches}+1):
            job = jd.Job(name=f"batch_{{i}}", function=lambda i=i: run_gromacs(i), cores={self.single_calculation_threads})
            jobs.append(job)

        dispatcher = jd.JobDispatcher(jobs, maxcores={self.total_threads}, engine="multiprocessing")
        results = dispatcher.run()

        # # Writing the outputs to a file
        # with open("orca_outputs.txt", "w") as outfile:
        #     for result in results.values():
        #         outfile.write(result)
        """
        )

    def postprocess(self):
        energies = []
        for batch_number in range(1, self.number_of_batches + 1):
            batchfile_path = os.path.join(self.calculator_folder, f"batch_{batch_number}")
            batch_energies = extract_energies(os.path.join(batchfile_path, "md.log"))
            print(batch_energies)
        # batch_number = 1
        # for i in range(1, self.settings.general.number_of_structures):
        #     threshold = batch_number * (self.structures_per_batch / self.settings.general.number_of_structures)
        #     if (i/self.settings.general.number_of_structures) < threshold:


def build_exports_string(input_string):
    """
    Constructs a bash export string from key-value pairs in the input string, to
    be inserted in the PBS/SLURM launch command.

    The input string should contain key-value pairs separated by semicolons (;),
    with each key and value separated by a colon (:). This function processes
    each key-value pair to create a string of bash export commands, appending
    the value of each key to its existing value (if any) in the environment.
    Each export statement is concatenated with '&&' to allow for sequential
    execution in bash.

    Args:
        input_string (str): A string containing key-value pairs in the format
                            "KEY1:VALUE1;KEY2:VALUE2", where each pair represents
                            an environment variable and its value to be exported.

    Returns:
        str: A string containing bash export commands for each key-value pair,
             suitable for execution in a shell script. Returns an empty string
             if the input is None or empty.

    Examples:
        >>> build_exports_string("PATH:/usr/local/bin;LD_LIBRARY_PATH:/usr/local/lib")
        'export PATH=/usr/local/bin:$PATH && export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH && '
    """
    # Check for None or empty input
    if not input_string:
        return ""

    # Split the input string into pairs
    pairs = input_string.split(";")

    # Initialize an empty list for export statements
    exports_string = ""

    for pair in pairs:
        # Split each pair by ':' and strip white spaces
        if ":" in pair:
            key, value = pair.split(":", 1)
            key = key.strip()
            value = value.strip()
            exports_string += f"export {key}={value}:${key} && "

    return exports_string


def build_modules_string(input_string):
    """
    Generates a string to load modules, formatted for bash, from a semicolon-separated list.

    Parses an input string containing module names separated by semicolons (;) and constructs
    a bash command string that loads these modules using `module load`. Each module load command
    is concatenated with '&&' for sequential loading.

    Args:
        input_string (str): A semicolon-separated list of module names to be loaded.

    Returns:
        str: A bash command string that loads each specified module. Returns an empty string
             if the input is None or empty.

    Example:
        >>> build_modules_string("gcc;python3.8;openmpi")
        'module load gcc && module load python3.8 && module load openmpi && '
    """
    # Check for None or empty input
    if not input_string:
        return ""

    # Split the string into pairs
    modules = input_string.split(";")

    modules_string = ""

    for module in modules:
        modules_string += f"module load {module.strip()} && "

    return modules_string
