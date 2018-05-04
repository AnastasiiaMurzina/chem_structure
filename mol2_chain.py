from penta_with_rotate import get_penta_points, find_section, find_section_and_rotate,\
    show_points, rotate_by_basis, show_named_points, get_reversed_section_and_basis, \
    show_named_points1, find_basis, find_basis_mave
from mol2_worker import xyz_names, xyz_names_bonds, Atom, atoms_and_bonds, Bond
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
eps_length = 0.001
pp = get_penta_points()


###########################Builder structure if possile else return fail########################
def prepare_bonds(bonds):
    for i, j in enumerate(bonds):
        if j[0] > j[1]:
            bonds[i][0], bonds[i][1] = bonds[i][1], bonds[i][0]
    return bonds


def bonds_to_one_way_dict(bonds):
    d = {}
    for i in bonds:
        if d.get(i[0]):
            d[i[0]].append(i[1::])
        else:
            d.update({i[0]: [i[1::]]})
    return d


def to_two_ways_bond2(one_way_bonds, with_attr=False):
    two_ways = {}
    if with_attr:
        for key, item in one_way_bonds.items():
            for i in sorted(item):
                if two_ways.get(key):
                    two_ways[key].append(i)
                else:
                    two_ways.update({key: [i]})
                if two_ways.get(i[0]):
                    two_ways[i[0]].append([key, i[1]])
                else:
                    two_ways.update({i[0]: [[key, i[1]]]})
        return two_ways
    for i in one_way_bonds:
        if two_ways.get(i[0]):
            two_ways.get(i[0]).append(i[1])
        else:
            two_ways.update({i[0]: [i[1]]})
        if two_ways.get(i[1]):
            two_ways.get(i[1]).append(i[0])
        else:
            two_ways.update({i[1]: [i[0]]})
    return two_ways


def coordinates_to_notation(info_from_file, valid_length=eps_length, save_length=False):
    '''
    :param info_from_file: read_from_file tuple
    :return: chain in dictionary notation (see above)
    '''
    bonds, atoms = info_from_file
    positions_copy = copy.deepcopy(atoms)
    notation = {}
    for key, item in atoms.items():
        cur_p = positions_copy.pop(key).position()
        connected = [k for k, ite in atoms.items() if abs(np.linalg.norm(ite.position() - cur_p) - 1) < valid_length]
        sections = []
        basis = find_basis(cur_p, [atoms[i].position() for i in connected])
        for i in connected:
            sections.append([i, find_section(cur_p, atoms[i].position(), basis0=basis, all_posibility=True)])
        notation.update({key: [sections, basis]})
    if save_length:
        bonds_keeper = {}
        bonds = bonds_to_one_way_dict(prepare_bonds(bonds))
        for key, item in notation.items():
            for i in item[0]:
                if key < i[0]:
                    attr = ''
                    for j in bonds[key]:
                        if j[0] == i[0]:
                            attr = j[1]
                            break
                    bond = Bond(key, i[0], attr, length=atoms[key].position()-atoms[i[0]].position(), section=i[1])
                    bonds_keeper.update({(key, i[0]): bond})
        return notation, bonds_keeper
    return notation


def mol2_to_notation(info_from_file):
    '''
    :param info_from_file: read_from_file tuple
    WARNING: may be keep Atoms without coordinates and Bonds with sections
    :return: chain in dictionary notation (see above)
    '''
    bonds, atoms = info_from_file
    positions_copy = copy.deepcopy(atoms)
    notation = {}
    bonds = bonds_to_one_way_dict(prepare_bonds(bonds))
    bonds2 = to_two_ways_bond2(bonds, with_attr=True)
    for key, item in atoms.items():
        cur_p = positions_copy.pop(key).position()
        connected = [i[0] for i in bonds2[key]]
        basis = find_basis(cur_p, [atoms[i].position() for i in connected])
        atoms[key].set_orientation(basis)
        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position(), basis0=basis, let_accurance=0.4)]
                                     for i in connected]), basis]})
    for key, item in bonds.items():
        for i in range(len(item)):
            bonds[key][i].insert(1, np.linalg.norm(atoms[key].position()-atoms[item[i][0]].position()))
    return notation, bonds


def write_mol_file(file_name, atoms, positions, bonds=[], attrs = {}):
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
        f1.write(str(len(atoms))+'\n\n')
        f1.write('@<TRIPOS>ATOM\n')
        for num, key in atoms.items():
            f1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(num+1, key.name, str(key.x),
                                                   str(key.y), str(key.z), key.name_i,
                                                                       key.i1, key.i2, key.i3))
        f1.write('@<TRIPOS>BOND\n')
        if bonds != []:
            for num, i in enumerate(bonds):
                # print(i, 't',attrs.get(tuple([i[0], i[1]])), attrs.get(tuple([i[1], i[0]])))
                f1.write("{0}\t{1}\t{2}\t{3}\n".format(str(num+1), str(i[0]),
                                                       str(i[1]), attrs[(i[0], i[1])][1] if attrs.get((i[0], i[1])) else attrs[(i[1], i[0])][1]))#str(i[2]) if len(i) > 2 else '1'))


def write_mol2_file(file_name, atoms, positions, bonds=[], attrs = {}):
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
        f1.write(str(len(atoms))+'\n\n')
        f1.write('@<TRIPOS>ATOM\n')
        for num, key in atoms.items():
            f1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(num+1, key.name, str(key.x),
                                                   str(key.y), str(key.z), key.name_i,
                                                                       key.i1, key.i2, key.i3))
        f1.write('@<TRIPOS>BOND\n')
        if bonds != []:
            for num, i in enumerate(bonds):
                # print(i, 't',attrs.get(tuple([i[0], i[1]])), attrs.get(tuple([i[1], i[0]])))
                f1.write("{0}\t{1}\t{2}\t{3}\n".format(str(num+1), str(i[0]),
                                                       str(i[1]), attrs[(i[0], i[1])][1] if attrs.get((i[0], i[1])) else attrs[(i[1], i[0])][1]))#str(i[2]) if len(i) > 2 else '1'))



##############################Mol2_preparations###########################################

def to_one_way_bonds(two_way_bonds):
    one_way = []
    for key, item in two_way_bonds.items():
        ix = 0
        if isinstance(item[0], (int, float)):
            ix = 1
        for i in item[ix]:
            if [i[0], key] not in one_way:
                one_way.append([key, i[0]])
    return one_way


def to_two_ways_bond(one_way_bonds, with_attr=False):
    two_ways = {}
    print(sorted(one_way_bonds))
    if with_attr:
        for i in one_way_bonds:
            if two_ways.get(i[0]):
                two_ways.get(i[0]).append((i[1], i[2]))
            else:
                two_ways.update({i[0]: [(i[1], i[2])]})
            if two_ways.get(i[1]):
                two_ways.get(i[1]).append((i[0], i[2]))
            else:
                two_ways.update({i[1]: [(i[0], i[2])]})
        return two_ways

    for i in one_way_bonds:
        if two_ways.get(i[0]):
            two_ways.get(i[0]).append(i[1])
        else:
            two_ways.update({i[0]: [i[1]]})
        if two_ways.get(i[1]):
            two_ways.get(i[1]).append(i[0])
        else:
            two_ways.update({i[1]: [i[0]]})
    return two_ways



if __name__ == '__main__':
######################mol2_files##############
    '''Test: read from mol2-file by atoms_and_bonds() - function
    allows us get information about every atom;
    xyz_names_bonds()- function
    '''
    atoms_info = atoms_and_bonds('Aniline.mol2')
    ln = mol2_to_notation(xyz_names_bonds('Aniline.mol2'))
    # ln = coordinates_to_notation(xyz_names_bonds('Aniline.mol2'), valid_length=0.5, save_length=True)
    # d31 = dimensional_length_unique_basis(ln)
    # for key, item in d31.items():
    #     atoms_info[key].set_xyz(item[0], item[1], item[2])
    # write_mol_file('My_aniline.mol2', atoms_info, ln)
    # ln = coordinates_to_mm_basis_notation(xyz_names_bonds('Caffein.mol2'), valid_length=0.5, save_length=True)

    # ln = coordinates_to_notation(xyz_names('benzene.mol2'), save_length=True)
    # write_mol_file('My_aniline.mol2', nms, d31, bonds=to_one_way_bonds(ln[0]), attrs=ln[2])
