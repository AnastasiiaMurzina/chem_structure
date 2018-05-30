import numpy as np
from penta_with_rotate import Ry, Rz
from mol2_worker import xyz_names_bonds
from mol2_chain import write_mol2_file


def get_rotated(file_name):
    return


if __name__ == '__main__':
    bonds, atoms = xyz_names_bonds('Caffein.mol2')
    # print(atoms)
    positions = {}
    for key, atom in atoms.items():
        positions.update({key: atom.position()})
    write_mol2_file('Rot_caffein.mol2', atoms, positions, bonds)
