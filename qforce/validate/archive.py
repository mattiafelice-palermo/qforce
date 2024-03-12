from openbabel import openbabel
import os
import shutil
from typing import List


def split_pdb_to_xyz(pdb_file_path, destination_folder_path):
    # Set the warning level to Error only to suppress warnings
    openbabel.obErrorLog.SetOutputLevel(openbabel.obError)

    # Open the PDB file
    if not os.path.isfile(pdb_file_path):
        print(f"The file {pdb_file_path} does not exist.")
        raise FileNotFoundError(f"Structures file does not exist: {pdb_file_path}")

    if os.path.exists(destination_folder_path):
        shutil.rmtree(destination_folder_path)
    os.makedirs(destination_folder_path)

    # Create a new OBMol object to handle the molecule
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "xyz")

    mol = openbabel.OBMol()
    if obConversion.ReadFile(mol, pdb_file_path):
        i = 1
        while mol:
            # Write the molecule to an XYZ file
            basename = f"molecule_{i}.xyz"
            obConversion.WriteFile(mol, os.path.join(destination_folder_path, basename))
            mol = openbabel.OBMol()  # Create a new OBMol to read the next molecule
            obConversion.Read(mol)
            i += 1
            if mol.NumAtoms() == 0:
                break
    else:
        print("Failed to read the PDB file.")


def count_xyz_files_in_folder(folder_path: str) -> int:
    """
    Counts the number of '.xyz' files in the specified folder.

    Args:
        folder_path (str): Path to the folder containing files.

    Returns:
        int: Number of '.xyz' files in the folder.
    """
    all_files: List[str] = os.listdir(folder_path)

    # Filter the list to include only files ending with ".xyz"
    xyz_files: List[str] = [
        name for name in all_files if os.path.isfile(os.path.join(folder_path, name)) and name.endswith(".xyz")
    ]

    return len(xyz_files)


def extract_energies(md_log_path):
    """
    Extracts energy information from a GROMACS md.log file and writes the formatted data to an output file.

    Parameters:
    md_log_path (str): The path to the md.log file.
    output_path (str): The path to the output file where the formatted energy information will be written.
    """
    # Flag to start capturing energy data
    capture = False
    frame_energies = []
    energies = []

    with open(md_log_path, "r") as log_file:
        for line in log_file:
            # Check if the current line indicates the start of energy data
            if line.strip().startswith("Energies (kJ/mol)"):
                capture = True
                continue

            # Check if capturing has ended based on a blank line following energy data
            if capture and line.strip() == "":
                capture = False
                energies.append(frame_energies)
                frame_energies = []
                continue

            if line.strip().startswith("Bond") or line.strip().startswith("Potential"):
                continue

            # If currently capturing energy data, process and store it
            if capture:
                splitted = line.strip().split()
                frame_energies.extend(splitted)
    return energies
