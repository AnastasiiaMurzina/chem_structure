import numpy as np
from penta_with_rotate import Ry, Rz
from mol2_worker import xyz_names_bonds
from mol2_chain import write_mol2_file, mol2_to_notation, bonds_of_paired


def get_rotated(name):
    # name = 'solution_neural_Ru_eaxial'
    bonds, atoms = xyz_names_bonds(name+'.mol2')
    ln = mol2_to_notation([bonds, atoms])
    positions = {}
    for key, atom in atoms.items():
        positions.update({key: atom.position().dot(Rz(2*np.pi/5))})
    write_mol2_file('Rot_'+name+'.mol2', atoms, positions, bonds_of_paired(ln[1]))


if __name__ == '__main__':
    name = 'solution_neural_Ru_eaxial'
    bonds, atoms = xyz_names_bonds(name+'.mol2')
    ln = mol2_to_notation([bonds, atoms])
    positions = {}
    for key, atom in atoms.items():
        positions.update({key: atom.position().dot(Rz(2*np.pi/5))})
    write_mol2_file('Rot_'+name+'.mol2', atoms, positions, bonds_of_paired(ln[1]))
