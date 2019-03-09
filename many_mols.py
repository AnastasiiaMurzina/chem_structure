import numpy as np
import copy
from mol2_worker import xyz_names_bonds
from mol2_chain_q import mol2_to_notation, to_two_ways_bond2, dimensional_structure, bonds_of_paired, to_two_ways_bond
from penta_with_rotate import find_section, rotate_by_basis, rotate_non_perpendicular, n_z, n_y, rotate_ten_vars, pp, find_basis


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

def large_order(not_atoms):
    '''
    :param not_atoms:
    :return: [(count_of_atoms, num_of_mol),]
    '''
    return sorted([(len(i) if not isinstance(i, int) else 1, k) for k, i in not_atoms.items()])


def get_notation_many_mols(atoms, bonds, method='first', **kwargs):
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
            ln = mol2_to_notation([lig_bs, lig_as], method=method, **kwargs)
            not_atoms.update({i: ln[0]})
            not_bonds.update({i: ln[1]})
        else:
            not_bonds.update({i: []})
            not_atoms.update({i: [ixx for ixx in lig_as.keys()][0]})
    return not_atoms, not_bonds


def places_for_zero_bonds(atoms, bonds, method='first', firmly=False, **kwargs):
    '''
    :param atoms:
    :param bonds: see get_notation_many_mols
    :return: {mol_to_connect: [[length, atom_with_connect, atom_of_this_mol],[...]], ...}
    '''
    not_atoms, not_bonds = get_notation_many_mols(atoms, bonds, method=method, **kwargs)
    # div_atoms, _, _ = molecular_divider(atoms, bonds)
    finder_zeros = large_order(not_atoms)
    nearests = {}
    for n_j, j in enumerate(finder_zeros[:-1:]):
        dss = []
        checked_atoms = not_atoms[j[1]]
        for _, i in finder_zeros[n_j+1:]:
            compared_atoms = not_atoms[i]
            if not isinstance(compared_atoms, int):
                if isinstance(checked_atoms, int):
                    checked_atoms = {checked_atoms: []}
                if firmly and len(checked_atoms.keys()) > 2:
                    checked_non_H_atoms = list(filter(lambda x: len(checked_atoms[x][0]) > 1, checked_atoms.keys()))
                else:
                    checked_non_H_atoms = checked_atoms.keys()
                for checked in checked_non_H_atoms:
                    cur_pos = atoms[checked].position()
                    for comp_atom in compared_atoms:
                        if not firmly or len(not_atoms[i][comp_atom][0]) > 1:
                            ds = np.linalg.norm(atoms[comp_atom].position() - cur_pos)
                            dss.append([ds, comp_atom, checked])
        nearests.update({j[1]: sorted(dss)[:2]})
    return nearests


def zero_connecter(atoms, to_zero_bonds, method='first', **kwargs):
    '''
    :param atoms: atoms from xyz_names_bonds(...)[1]
    :param to_zero_bonds: places_for_zero_bonds
    :return: {mol: [[length, this, another_mol], section], ...}
    '''
    zero_bonds = {}
    for key, item in to_zero_bonds.items():
        zero_bonds.update({key: []})
        for i in item:
            zero_bonds[key].append([i, find_section(atoms[i[1]].position(), atoms[i[2]].position(), basis0=find_basis(atoms[i[1]].position(), atoms[i[2]].position(), method=method, kwargs=kwargs), method=method, kwargs=kwargs)])
    return zero_bonds


def unpack_with_zero_bonds(atoms_not, bonds_not, zero_bonds, method='first', **kwargs):
    the_largest_flag = True
    def av():
        positions, avs = np.zeros(3), np.zeros(3)
        for i, sec in zero_bonds[key]:
            positions += (pp[sec] * i[0] + dim_structure[i[1]])
            avs += ds.get(i[2], np.zeros(3))
        return positions/2-avs/2
    for _, key in reversed(large_order(atoms_not)):
        if bonds_not[key] != []:
            ds = dimensional_structure([atoms_not[key], bonds_of_paired(bonds_not[key])], method=method, **kwargs)
            if the_largest_flag:
                the_largest_flag = False
                dim_structure = ds
            else: # if it isn't the largest one
                mover_vector = av()
                for k, i in ds.items():
                    dim_structure.update({k: i+mover_vector})
        else:
            dim_structure.update({atoms_not[key]: av()})
    return dim_structure


def write_mol2_file(file_name, atoms, positions, bonds):
    '''
    :param file_name:
    :param names: chemical elements names
    :param positions: xyz - coordinates
    :param bonds: one way bonds
    :return: None, void write function
    '''
    with open(file_name, 'w') as f1:
        f1.write('@<TRIPOS>MOLECULE\n')
        f1.write('some info\n')
        f1.write(str(len(atoms))+'\t'+str(len(bonds)-1)+'\n\n')
        f1.write('@<TRIPOS>ATOM\n')
        for num, key in atoms.items():
            f1.write("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(num, key.name, str(positions[num][0]),
                                                                            str(positions[num][1]), str(positions[num][2]),
                                                                            key.name_i, key.i1, key.i2, key.i3))
        f1.write('@<TRIPOS>BOND\n')
        for k, num in enumerate(bonds.items()):
            num, i = num
            for ix in i:
                if num < ix[0]:
                    # print(str(k+1), str(num), str(ix[0]), str(ix[1]))
                    f1.write("\t{0}\t{1}\t{2}\t{3}\n".format(str(k+1), str(num), str(ix[0]), str(ix[1])))


if __name__ == '__main__':

    # name = 'vacuum_cation_singlet_Fe_full'
    name = '/home/anastasiia/PycharmProjects/chem_structure/TS/TS-2c-2d_last'
    # name = '/home/anastasiia/PycharmProjects/chem_structure/TS/2b-2c/TS-2b-2c_step_16'
    bs, ass = xyz_names_bonds(name + '.mol2')
    div_atoms, ligs_as, lig_bonds = molecular_divider(ass, bs)
    atoms_notation, bonds_notation = get_notation_many_mols(ass, bs, method='ten')
    needs_to_zero_discribe = places_for_zero_bonds(ass, bs, firmly=False, method='ten')
    zero_bonds = zero_connecter(ass, needs_to_zero_discribe, method='ten')
    coords = unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds, method='ten')
    # print(needs_to_zero_discribe)
    # print(coords)
    write_mol2_file('2c-2d_my.mol2', ass, coords, to_two_ways_bond(bs, with_attr=True))
    ##################################################
    # import sys, os, subprocess, shutil, tempfile
    # d = os.getcwd()
    # tmpdir = tempfile.mkdtemp()
    #
    # compared_molecule2 = os.path.join(d, name + '.mol2')
    # compared_molecule = os.path.join(d, 'long.mol2')
    #
    # sim_file = os.path.join(tmpdir, 'sim_three.txt')
    # subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule2, compared_molecule, sim_file],
    #                 stdout=subprocess.DEVNULL)
    # with open(sim_file, 'r') as f:
    #     f.readline()
    #     coeffs = f.readline().split()
    #     coeffs = coeffs[2:6] #+ coeffs[8:]
    #     # print(coeffs)
    #     print(np.linalg.norm(np.array(list(map(float, coeffs)))))
    # shutil.rmtree(tmpdir)

    # print(unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds))
    # mol_with_length = sorted([(len(i) if not isinstance(i, int) else 1, k) for k, i in atoms_notation.items()])

    # for i in mol_with_length[:-1:]:
    #     print(i)

