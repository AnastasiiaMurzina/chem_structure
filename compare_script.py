import sys, os, subprocess, shutil

from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file

if __name__ == '__main__':
    compared_molecul = 'Aniline'
    atoms_info = atoms_and_bonds(compared_molecul + '.mol2')
    ln = mol2_to_notation(xyz_names_bonds(compared_molecul + '.mol2'))
    paired = bonds_of_paired(ln[1])
    dim_structure = dimensional_structure([ln[0], paired])
    write_mol2_file('My_' + compared_molecul + '.mol2', atoms_info, dim_structure, bonds=paired)