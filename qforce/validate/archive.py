from openbabel import openbabel
import os
import shutil
import time


def split_pdb_to_xyz(pdb_file_path, destination_folder_path):
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
