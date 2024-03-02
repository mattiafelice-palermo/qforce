from colt import Colt
from abc import ABC, abstractmethod
import shutil
import os
import textwrap
import subprocess

from pprint import pprint


def get_calculator(settings):
    """
    Returns the appropriate calculator based on the specified settings.

    Args:
        settings: The settings object containing the configurations for qforce-validate

    Returns:
        Calculator: An instance of a calculator based on the specified options.

    Raises:
        NotImplementedError: If the calculator specified in the settings is not supported.
    """
    scheduler = settings.general.calculator.lower()

    if scheduler == "orca":
        return Orca(settings)
    else:
        raise NotImplementedError(f"Calculator '{scheduler}' is not implemented.")


class CalculatorABC(ABC):
    def __init__(self, settings):
        self.settings = settings
        self._check_files_in_pool()

    def _check_files_in_pool(self):
        print("here")
        pool_path = self.settings.general.pool_path
        for file in os.listdir(pool_path):
            if file.endswith(".xyz"):
                return

        raise FileNotFoundError("No xyz file has been found in the pool directory.")

    @abstractmethod
    def run(self): ...


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

    def __init__(self, settings):
        super().__init__(settings)
        self.charge = settings.orca.charge
        self.multiplicity = settings.orca.multiplicity
        self.method = settings.orca.method
        self.blocks = settings.orca.blocks
        self.total_threads = settings.orca.total_threads
        self.total_memory = settings.orca.total_memory
        self.single_calculation_threads = settings.orca.single_calculation_threads
        self.orca_executable = settings.orca.orca_executable
        self.calculator_folder = os.path.join(self.settings.general.job_dir, self.settings.general.calculator)
        self.modules_string = build_modules_string(settings.orca.modules)
        self.exports_string = build_exports_string(settings.orca.exports)
        self.launch_command = "python dispatcher.py"
        self.queue = settings.orca.queue
        self.conda_environment = settings.orca.conda_environment
        self.job_dir = self.calculator_folder  # just for the schedulers

    def run(self, dry_run=False):
        self._setup_working_folder()
        self._generate_orca_template()
        self._generate_launch_script()

        if not dry_run:
            print("Running in system!")
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
        # Fill in the template with specific details
        # This is where the script is customized based on the user input
        number_of_structures = self.settings.general.number_of_structures + 1

        return textwrap.dedent(
            f"""
        import jobdispatcher as jd
        import subprocess
        import os
        import shutil
        
        def run_orca(filename):
            # Call the orca_calc.sh script and capture its output
            folder_name = create_directory(os.path.join("{self.calculator_folder}", filename))
            shutil.copy2('orca_template.inp', os.path.join(folder_name, f"{{filename}}.inp"))
            xyz_path = os.path.join("{self.settings.general.pool_path}", f"{{filename}}.xyz")
            subprocess.run(f"sed -i 's#FILENAME#{{xyz_path}}#g' {{filename}}.inp", shell=True, cwd=folder_name, check=True)

            result = subprocess.run(f'{self.modules_string}{self.exports_string} $(which orca) {{filename}}.inp > {{filename}}.out 2> {{filename}}.err', shell=True, cwd=folder_name, check=True)

            print(result)
            
            return result.stdout

        def create_directory(folder):
            if os.path.exists(folder):
                shutil.rmtree(folder)
            os.makedirs(folder)
            return folder
            
        jobs = []
        for i in range(1, {number_of_structures}):
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

        template_string = textwrap.dedent(
            f"""
        ! {self.method} PAL{self.single_calculation_threads} 
        * xyzfile {self.charge} {self.multiplicity} FILENAME
        """
        )

        with open(orca_input_path, "w") as orca_input_handle:
            orca_input_handle.write(template_string)


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
