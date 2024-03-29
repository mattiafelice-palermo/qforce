import os
import logging
from datetime import datetime
import traceback
import yaml
import textwrap

# =======================================================================================
# MISC.PY - Miscellaneous Utilities
# =======================================================================================
# This section of the misc.py module is dedicated to miscellaneous utility classes and
# functions that don't necessarily fit into the main thematic modules of the application
# but are essential for handling specific tasks such as file manipulation, string
# processing, or any other utility operations needed across the project.
#
# Included in this section:
# - GroBoxEditor: A class for editing and manipulating box vectors in GRO files. This
#   utility allows for reading, modifying, and retrieving box vector information from
#   GRO files, facilitating the manipulation of molecular simulation box dimensions.
#
# =======================================================================================
# Begin Miscellaneous Utilities
# =======================================================================================


class GroBoxEditor:
    """
    A class for editing the box vectors found in the last line of a GRO file.

    Attributes:
    file_path (str): Path to the GRO file.
    original_box_vectors (list): The original box vectors as extracted from the file.
    box_vectors (list): The box vectors that may be modified by the user.
    decimal_places (list): The number of decimal places for each box vector component.
    """

    def __init__(self, file_path):
        """Loads box vectors from a GRO file."""

        self.file_path = file_path
        self.original_box_vectors, self.decimal_places = self._load_box_vectors()
        self.box_vectors = self.original_box_vectors.copy()

    def _load_box_vectors(self):
        """Returns box vectors and their decimal places from the file's last line."""

        with open(self.file_path, "r") as file:
            lines = [line for line in file.readlines() if line.strip()]
            last_line = lines[-1].strip()
            vectors_str = last_line.split()
            vectors = [float(value) for value in vectors_str]
            decimal_places = [len(value.split(".")[1]) if "." in value else 0 for value in vectors_str]
            return vectors, decimal_places

    def edit_box_vector(self, x=None, y=None, z=None):
        """Edits x, y, z components of the box vector. Unspecified components remain unchanged."""

        if x is not None:
            self.box_vectors[0] = x
        if y is not None:
            self.box_vectors[1] = y
        if z is not None:
            self.box_vectors[2] = z

    def _format_vector_string(self, vectors):
        """Formats vectors into a string, maintaining original decimal precision."""

        return "   ".join(f"{vectors[i]:.{self.decimal_places[i]}f}" for i in range(3))

    def get_original_box_vector_string(self):
        """Returns modified box vectors as a string with preserved decimal precision."""
        return self._format_vector_string(self.original_box_vectors)

    def get_box_vector_string(self):
        """
        Retrieves the (potentially modified) box vectors as a formatted string, preserving the original decimal places.

        Returns:
        str: The current box vectors formatted as a string.
        """
        return self._format_vector_string(self.box_vectors)


# ==============================================================================
# PATH STRING MANIPULATION FUNCTIONS
# ==============================================================================
# This section contains functions dedicated to the manipulation and resolution
# of file and directory paths. These utilities include, but are not limited to,
# removing quotes from paths, expanding environment variables, resolving relative
# paths to absolute paths, and ensuring path validity and accessibility. These
# functions are designed to handle a variety of path specifications, making
# them versatile for use across different modules and applications where
# path string manipulation is required.
# ==============================================================================

# [Place your path manipulation functions here, like `remove_quotes` and `get_fullpath`]


def get_fullpath(raw_filepath, folder_path=None):
    """
    Resolves the full path of a structure file based on different path specifications,
    handling paths enclosed in quotes, expanding environment variables, resolving symbolic links,
    and checking for file existence and read permissions.

    Parameters:
    - filepath: The file path as specified, possibly with surrounding quotes.

    Returns:
    - The absolute path to the structure file, or None if the file does not exist or is not accessible.
    """
    # Remove surrounding quotes and expand environment variables from the file path
    filepath = remove_quotes(raw_filepath)
    filepath = os.path.expandvars(filepath)

    # Expand the user's home directory if the path starts with ~
    filepath = os.path.expanduser(filepath)

    # Resolve symbolic links
    filepath = os.path.realpath(filepath)

    # Check if the path is already absolute, if not, make it absolute
    if folder_path is None:
        cwd_path = os.getcwd()
    else:
        cwd_path = get_fullpath(folder_path)  # a bit recursive, but should work

    if not os.path.isabs(filepath):
        filepath = os.path.join(cwd_path, filepath)
    else:
        filepath = os.path.abspath(filepath)

    # Check for file existence and read permissions
    if not os.path.exists(filepath):
        print(f"Error: The file '{filepath}' does not exist.")
        return None
    if not os.access(filepath, os.R_OK):
        print(f"Error: The file '{filepath}' is not readable.")
        return None

    # Return the resolved absolute path
    return filepath


def remove_quotes(input_str):
    """
    Removes surrounding quotes from a string.

    Supports single ('), double ("), and backtick (`) quotes.

    Parameters:
    - input_str: The input string from which to remove quotes.

    Returns:
    - The string without surrounding quotes.
    """
    quote_chars = "'\"`´"  # Include all types of quote characters you want to handle
    if len(input_str) >= 2 and input_str[0] in quote_chars and input_str[-1] in quote_chars:
        return input_str[1:-1]
    return input_str


# ==============================================================================
# END OF PATH STRING MANIPULATION FUNCTIONS
# ==============================================================================

# ==============================================================================
# LOGGING FORMATTER
# ==============================================================================
# This section contains the class definition for the YAML formatter for the
# logging facility.
# ==============================================================================


class EnhancedYamlFormatter(logging.Formatter):
    last_timestamp = ""
    timestamp_counter = 0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        yaml.add_representer(str, self.representer_multiline_str, Dumper=yaml.SafeDumper)

    @staticmethod
    def representer_multiline_str(dumper, data):
        if "\n" in data:
            return dumper.represent_scalar("tag:yaml.org,2002:str", data, style="|")
        return dumper.represent_scalar("tag:yaml.org,2002:str", data)

    def format(self, record):
        timestamp_format = "%Y-%m-%d %H:%M:%S"
        formatted_timestamp = datetime.utcfromtimestamp(record.created).strftime(timestamp_format)

        # Check if this timestamp is the same as the last one
        if formatted_timestamp == self.last_timestamp:
            self.timestamp_counter += 1
        else:
            self.last_timestamp = formatted_timestamp
            self.timestamp_counter = 1  # Reset counter for new timestamp

        # Append the counter to the timestamp to make it unique
        unique_timestamp = f"{formatted_timestamp}.{self.timestamp_counter}"

        structured_log = {
            "time": unique_timestamp,
            "level": record.levelname,
            "line": record.lineno,
            "function": record.funcName,
            "file": record.filename,
            "message": record.getMessage(),
        }

        if record.exc_info:
            tb_list = traceback.format_exception(*record.exc_info)
            indented_tb = textwrap.indent("".join(tb_list), prefix="    ")
            structured_log["exception"] = indented_tb

        log_entry = {unique_timestamp: structured_log}

        return yaml.dump(log_entry, Dumper=yaml.SafeDumper, default_flow_style=False, sort_keys=False)


# class EnhancedYamlFormatter(logging.Formatter):
#     def format(self, record):
#         # Create a structured log record
#         structured_log = {
#             "time": datetime.utcfromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S"),
#             "level": record.levelname,
#             "line": record.lineno,
#             "function": record.funcName,
#             "file": record.filename,
#             "message": record.getMessage(),  # Gets the log message
#         }

#         # Check if exception information is included
#         if record.exc_info:
#             # Extract traceback as a string
#             tb_str = traceback.format_exception(*record.exc_info)
#             # Indent each line for YAML formatting
#             indented_tb = textwrap.indent("".join(tb_str), "    ")
#             structured_log["exception"] = "".join(indented_tb)

#         # Create a parent key based on the log creation time to maintain hierarchy
#         parent_key = datetime.utcfromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
#         log_entry = {parent_key: structured_log}

#         # Return the YAML representation of the log entry
#         return yaml.dump(log_entry, default_flow_style=False, sort_keys=False, allow_unicode=True)


# class EnhancedYamlFormatter(logging.Formatter):
#     def format(self, record):
#         # Create a dictionary with only the desired information
#         structured_log = {
#             "time": datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S"),
#             "level": record.levelname,
#             "line": record.lineno,
#             "function": record.funcName,
#             "file": record.filename,
#             "message": record.msg,
#         }

#         # Create a parent key based on the log creation time to maintain hierarchy
#         parent_key = datetime.fromtimestamp(record.created).strftime("%Y-%m-%d %H:%M:%S")
#         log_entry = {parent_key: structured_log}

#         # Return the YAML representation of the log entry
#         return yaml.dump(log_entry, default_flow_style=False, sort_keys=False)


# ==============================================================================
# END OF LOGGING FORMATTER CLASS DEFINITION
# ==============================================================================
