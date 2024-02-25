from openbabel import pybel
import os
import shutil


def split_pdb_to_xyz(pdb_file_path, destination_folder_path):
    # Open the PDB file
    molecules = pybel.readfile("pdb", pdb_file_path)

    if os.path.exists(destination_folder_path):
        shutil.rmtree(destination_folder_path)
    os.makedirs(destination_folder_path)

    # Iterate through each molecule and save as an XYZ file
    for i, mol in enumerate(molecules):
        basename = f"molecule_{i+1}.xyz"
        mol.write("xyz", os.path.join(destination_folder_path, basename), overwrite=True)
