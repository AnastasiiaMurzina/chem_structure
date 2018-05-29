import numpy as np
import copy
from mol2_worker import xyz_names_bonds
from mol2_chain import mol2_to_notation, to_two_ways_bond2, dimensional_structure, bonds_of_paired
from penta_with_rotate import find_section, rotate_by_basis, rotate_non_perpendicular, n_z, n_y, rotate_ten_vars, pp


def molecular_divider(atoms, bonds):
    '''
    :param atoms: bonds, atoms = xyz_names_bonds()
    :return: {atom_num: mol_num} (prob. later: div_atoms - dictionary like {mol_num: [atoms_of_this_molecule], ...},)
     [lig_as] - list of dictionaries with atoms of structure
     [lig_bs] - list of dictionaries with bonds of structure
    '''
    atom_names = set([j.i2 for _, j in atoms.items()])
    div_atoms, ligs_atoms, ligs_bonds = {}, {}, {}
    for ix, atom in enumerate(atom_names):
        lig_as = dict(filter(lambda x: x[1].i2 == atom, atoms.items()))
        lig_bs = list(filter(lambda x: x[0] in lig_as.keys() or x[1] in lig_as.keys(), bonds))
        for lig in lig_as:
            div_atoms.update({lig: ix})
        ligs_atoms.update({ix: lig_as})
        ligs_bonds.update({ix: lig_bs})
    return div_atoms, ligs_atoms, ligs_bonds


def get_notation_many_mols(atoms, bonds):
    '''
    :param atoms:
    :param bonds: bonds, atoms = xyz_names_bonds()
    :return: dictionary of notations for every molecule
    for atoms and for bonds separately
    '''
    div_atoms, ligs_as, lig_bonds = molecular_divider(atoms, bonds)
    not_atoms, not_bonds = {}, {}
    for i in ligs_as.keys():
        lig_as = ligs_as[i]
        lig_bs = lig_bonds[i]
        if lig_bs != []:
            ln = mol2_to_notation([lig_bs, lig_as])
            not_atoms.update({i: ln[0]})
            not_bonds.update({i: ln[1]})
        else:
            not_bonds.update({i: []})
            not_atoms.update({i: [ixx for ixx in lig_as.keys()][0]})
    return not_atoms, not_bonds


def places_for_zero_bonds(atoms, bonds):
    '''
    :param atoms:
    :param bonds: see get_notation_many_mols
    :return: {mol_to_connect: [[length, atom_with_connect, atom_of_this_mol],[...]], ...}
    '''
    not_atoms, not_bonds = get_notation_many_mols(atoms, bonds)
    # div_atoms, _, _ = molecular_divider(atoms, bonds)
    finder_zeros = sorted([(len(i) if not isinstance(i, int) else 1, k) for k, i in not_atoms.items()])
    nearests = {}
    for n_j, j in enumerate(finder_zeros[:-1:]):
        dss = []
        checked_atoms = not_atoms[j[1]]
        for _, i in finder_zeros[n_j+1:]:
            compared_atoms = not_atoms[i]
            if not isinstance(compared_atoms, int):
                if isinstance(checked_atoms, int):
                    checked_atoms = {checked_atoms: []}
                if len(checked_atoms.keys()) > 2:
                    checked_non_H_atoms = list(filter(lambda x: len(checked_atoms[x][0]) > 1, checked_atoms.keys()))
                else:
                    checked_non_H_atoms = checked_atoms.keys()
                for checked in checked_non_H_atoms:
                    cur_pos = atoms[checked].position()
                    for comp_atom in compared_atoms:
                        if len(not_atoms[i][comp_atom][0]) > 1:
                            ds = np.linalg.norm(atoms[comp_atom].position() - cur_pos)
                            dss.append([ds, comp_atom, checked])
        nearests.update({j[1]: sorted(dss)[:2]})
    print(nearests)
    return nearests


def zero_connecter(atoms, to_zero_bonds):
    '''
    :param atoms: atoms from xyz_names_bonds(...)[1]
    :param to_zero_bonds: places_for_zero_bonds
    :return: {mol: [[length, this, another_mol], section], ...}
    '''
    zero_bonds = {}
    for key, item in to_zero_bonds.items():
        zero_bonds.update({key: []})
        for i in item:
            zero_bonds[key].append([i, find_section(atoms[i[1]].position(), atoms[i[2]].position())])
    return zero_bonds



def unpack_with_zero_bonds(atoms_not, bonds_not, zero_bonds):
    for key, item in atoms_not.items():
        if bonds_not[key] != []:
            ds = dimensional_structure([item, bonds_of_paired(bonds_not[key])])
            if zero_bonds.get(key):
                print(key, zero_bonds[key])
        # print(item)
        # print(bonds_not[key])
    dim_structure = []
    return dim_structure


if __name__ == '__main__':
    name = 'vacuum_cation_singlet_Fe_full'
    bs, ass = xyz_names_bonds(name + '.mol2')
    div_atoms, ligs_as, lig_bonds = molecular_divider(ass, bs)
    atoms_notation, bonds_notation = get_notation_many_mols(ass, bs)
    # print(atoms_notation)
    # print('boonds')
    # print(bonds_notation)
    needs_to_zero_discribe = places_for_zero_bonds(ass, bs)
    # print(needs_to_zero_discribe)
    zero_bonds=zero_connecter(ass, needs_to_zero_discribe)
    unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds)
    # mol_with_length = sorted([(len(i) if not isinstance(i, int) else 1, k) for k, i in atoms_notation.items()])

    # for i in mol_with_length[:-1:]:
    #     print(i)

