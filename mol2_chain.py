from penta_with_rotate import get_penta_points, find_section, find_section_and_rotate,\
    show_points, rotate_by_basis, show_named_points, get_reversed_section_and_basis, \
    show_named_points1, find_basis, find_basis_mave
from mol2_worker import xyz_names, xyz_names_bonds
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import copy

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
eps_length = 0.001
pp = get_penta_points()


###########################Builder structure if possile else return fail########################
def dimensional_structure(chain_v):
    '''
    :param chain_v: dictionary of atoms with bonds and atoms' basis properties
    :return: list of atoms with their center coordinates if it's possible else None
    '''
    dim_structure = {1: [np.array([0, 0, 0]), np.array([0, 0])]}
    chain_copy = copy.deepcopy(chain_v)
    p = chain_copy.pop(1)
    p.insert(0, 1) #p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        for i in bonds: # build bonds for cur_key atom
            if i[0] in dim_structure: # if we already have position: we'll check
                section, b0, b1 = find_section_and_rotate(dim_structure[cur_key][0], dim_structure[i[0]][0])
                ix = [j for j in chain_v[cur_key][0] if j[0]==i[0]][0]
                if (ix[1], chain_v[cur_key][1][0], chain_v[cur_key][1][1]) != (section, b0, b1):
                    print('Invalid structure, cause {0} is not {1} for atom {2}'.format(str(ix[1]),
                                                                                        str(section), str(i[0])))
                    return
            else: #if it's new bond
                coord = rotate_by_basis(pp[i[1]], dim_structure[cur_key][1][0], dim_structure[cur_key][1][1])+dim_structure[cur_key][0]
                section, angle0, angle1 = find_section_and_rotate(dim_structure[cur_key][0], coord)
                dim_structure.update({i[0]: [coord, [angle0, angle1]]})
                poper = chain_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure, chain_v


def saver_l(info_from_file, notation, bonds_attr):
    positions = info_from_file[0]
    lengths = {}
    for key, item in notation.items():
        for i in item[0]:
            ck = i[0] if isinstance(i[0], (float, int)) else i[0][0]
            if bonds_attr.get((key, ck)) or bonds_attr.get(((ck, key))):
                lengths.update({(key, ck): [np.linalg.norm(positions[key]-positions[int(ck)]), bonds_attr[(key, ck)] if bonds_attr.get((key, ck)) else bonds_attr[(ck, key)]]})
    return lengths


def coordinates_to_notation(info_from_file, valid_length=eps_length, save_length=False):
    '''
    :param info_from_file: read_from_file tuple
    :return: chain in dictionary notation (see above)
    '''
    positions, names = info_from_file
    positions_copy = copy.deepcopy(positions)
    notation = {}
    for key, item in positions.items():
        cur_p = positions_copy.pop(key)
        connected = [k for k, ite in positions.items() if abs(np.linalg.norm(ite - cur_p) - 1) < valid_length]
        sections = []
        basis = np.zeros(2)
        for i in connected:
            if key == 1:
                # sections.append((find_section(cur_p, positions[i]), 0, 0))
                sections.append((find_section_and_rotate(cur_p, positions[i]), 0, 0))
            else:
                s, b0, b1 = find_section_and_rotate(cur_p, positions[i])[::]
                basis = np.array([b0, b1])
                sections.append(find_section_and_rotate(cur_p, positions[i]))
        app = []
        for inx in range(len(connected)):
            app.append([connected[inx], sections[inx][0]])
        notation.update({key: [[a for a in app], basis]})
#########################################################
    upd = []
    if notation[1][0][0][1] == None:
        for num, i in enumerate(notation[1][0]):
            if i[1] == None:
                upd.append([i[0], get_reversed_section_and_basis(notation[i[0]][0][0][1], notation[i[0]][1])])
        # bss = []
        # for i in upd:
        #     bss.append(np.array([i[1][1], i[1][2]]))
        notation.update({1: [[[i[0], i[1]] for i in upd], np.zeros(2)]})
    if save_length: return notation, names, saver_l(info_from_file, notation)
##################################################
    return notation, names


def dimensional_structure_from_efa_notation(notation):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    :return: xyz-info
    '''
    dim_structure = {1: [np.array([0, 0, 0]), np.array([0, 0])]}
    bonds_l, name = notation
    p = bonds_l[1]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, 1)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        # print(cur_key, bonds, basis)
        for i in bonds:  # build bonds for cur_key atom
            if i[0] in dim_structure:  # if we already have position: we'll check
                if isinstance(bonds_l[i[0]][1][0], (float, int)):
                    section, b0, b1 = find_section_and_rotate(bonds_l[cur_key][0], bonds_l[i[0]][0])
                    ix = [j for j in bonds_l[cur_key][0] if j[0] == i[0]][0]
                    if (ix[1], bonds_l[cur_key][1][0], bonds_l[cur_key][1][1]) != (section, b0, b1):
                        print('Invalid structure, cause {0} is not {1} for atom {2}'.format(str(ix[1]),
                                                                                            str(section), str(i[0])))
                        return
            else:  # if it's new bond
                if isinstance(i[1], (float, int)):
                    coord = rotate_by_basis(pp[i[1]], dim_structure[cur_key][1][0], dim_structure[cur_key][1][1]) + \
                        dim_structure[cur_key][0]
                    section, angle0, angle1 = find_section_and_rotate(dim_structure[cur_key][0], coord)
                    dim_structure.update({i[0]: [coord, [angle0, angle1]]})
                else:
                    coord = rotate_by_basis(pp[i[1][0]], i[1][1], i[1][2]) + \
                            dim_structure[cur_key][0]
                    section, angle0, angle1 = find_section_and_rotate(dim_structure[cur_key][0], coord)
                    dim_structure.update({i[0]: [coord, [angle0, angle1]]})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure, name


def dimensional_from_lefa_notation(notation):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    with length
    :return: xyz-info
    '''
    dim_structure = {1: [np.array([0, 0, 0]), np.array([0, 0])]}
    bonds_l, name, lengths = notation
    p = bonds_l[1]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, 1)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        # print(cur_key, bonds, basis)
        for i in bonds:  # build bonds for cur_key atom
            if i[0] in dim_structure:  # if we already have position: we'll check
                if isinstance(bonds_l[i[0]][1][0], (float, int)):
                    section, b0, b1 = find_section_and_rotate(bonds_l[cur_key][0], bonds_l[i[0]][0])
                    ix = [j for j in bonds_l[cur_key][0] if j[0] == i[0]][0]
                    if (ix[1], bonds_l[cur_key][1][0], bonds_l[cur_key][1][1]) != (section, b0, b1):
                        print('Invalid structure, cause {0} is not {1} for atom {2}'.format(str(ix[1]),
                                                                                            str(section), str(i[0])))
                        return
            else:  # if it's new bond
                if isinstance(i[1], (float, int)):
                    coord = rotate_by_basis(pp[i[1]], dim_structure[cur_key][1][0], dim_structure[cur_key][1][1])*lengths[(cur_key, i[0])] + \
                        dim_structure[cur_key][0]
                    section, angle0, angle1 = find_section_and_rotate(dim_structure[cur_key][0], coord)
                    dim_structure.update({i[0]: [coord, [angle0, angle1]]})
                else:
                    coord = rotate_by_basis(pp[i[1][0]], i[1][1], i[1][2])*lengths[(cur_key, i[0])] + \
                            dim_structure[cur_key][0]
                    section, angle0, angle1 = find_section_and_rotate(dim_structure[cur_key][0], coord)
                    dim_structure.update({i[0]: [coord, [angle0, angle1]]})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure, name


def dimensional_length_unique_basis(notation):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    with length
    :return: xyz-info
    '''
    dim_structure = {1: np.array([0, 0, 0])}#, np.array([0, 0])]}
    bonds_l, name, lengths = notation
    p = bonds_l[1]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, 1)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        for i in bonds:  # build bonds for cur_key atom
            if not (i[0] in dim_structure):  # if we don't have position:
                coord = rotate_by_basis(pp[i[1]], basis[0], basis[1])*lengths[(cur_key, i[0])][0]+dim_structure[cur_key]
                dim_structure.update({i[0]: coord})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure, name


def get_bonds_attr(bonds):
    attr = {}
    for i in bonds:
        attr.update({(i[0], i[1]): i[2]})
    return attr

def coordinates_to_mm_basis_notation(info_from_file, valid_length=eps_length, save_length=False):
    '''
    :param info_from_file: read_from_file tuple
    :return: chain in dictionary notation (see above)
    '''
    positions, names, bonds = info_from_file
    positions_copy = copy.deepcopy(positions)
    notation = {}
    for key, item in positions.items():
        cur_p = positions_copy.pop(key)
        connected = [k for k, ite in positions.items() if abs(np.linalg.norm(ite - cur_p) - 1) < valid_length]
        sections = []
        basis = find_basis(cur_p, [positions[i] for i in connected])
        for i in connected:
            sections.append([i, find_section(cur_p, positions[i], basis0=basis, all_posibility=True)])
        notation.update({key: [sections, basis]})
    if save_length:
        return notation, names, saver_l(info_from_file, notation, get_bonds_attr(bonds))
    return notation, names


def coordinates_to_min_mean_basis_notation(info_from_file, valid_length=eps_length, save_length=False):
    '''
    Unique basis for one node by minimum of average distances
    :param info_from_file: read_from_file tuple
    :return: chain in dictionary notation (see above)
    '''
    positions, names = info_from_file
    positions_copy = copy.deepcopy(positions)
    notation = {}
    for key, item in positions.items():
        cur_p = positions_copy.pop(key)
        connected = [k for k, ite in positions.items() if abs(np.linalg.norm(ite - cur_p) - 1) < valid_length]
        sections = []
        basis = find_basis_mave(cur_p, connected)
        for i in connected:
            sections.append([i, find_section(cur_p, positions[i], basis0=basis, all_posibility=True)])
        notation.update({key: [sections, basis]})
    if save_length:
        return notation, names, saver_l(info_from_file, notation)
    return notation, names


def mol2_coordinates_to_notation(info_from_file, save_length=False):
    '''
    :param info_from_file: read_from_file tuple
    :return: chain in dictionary notation (see above)
    '''
    positions, names, bonds = info_from_file
    positions_copy = copy.deepcopy(positions)
    notation = {}
    for key, item in positions.items():
        cur_p = positions_copy.pop(key)
        connected = [k for k, ite in positions.items() if abs(np.linalg.norm(ite - cur_p) - 1) < 0.3]
        sections = []
        basis = find_basis(cur_p, [positions[i] for i in connected])
        for i in connected:
            sections.append([i, find_section(cur_p, positions[i], basis0=basis, all_posibility=True)])
        notation.update({key: [sections, basis]})
    if save_length:
        return notation, names, saver_l(info_from_file, notation)
    return notation, names

##############################################################################################


def write_xyz_file(file_name, names, positions):
    '''
    :param file_name:
    :param names: chemical elements names
    :param positions: xyz - coordinates
    :return: None, void write function
    '''
    with open(file_name, 'w') as f1:
        f1.write(str(len(names))+'\n\n')
        for key, item in names.items():
            f1.write("{0}\t{1}\t{2}\t{3}\n".format(item, str(positions[key][0]),str(positions[key][1]),str(positions[key][2])))


def read_xyz_file(file):
    with open(file, 'r') as file:
        n = int(file.readline())
        names = {}
        positions = {}
        file.readline()
        for i in range(n):
            line = file.readline().split()
            names.update({i+1: line[0]})
            x, y, z = [float(j) for j in line[1::]]
            positions.update({i+1: np.array([x,y,z])})
    return positions, names


def write_mol_file(file_name, names, positions, bonds=[], attrs = {}):
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
        f1.write(str(len(names))+'\n\n')
        f1.write('@<TRIPOS>ATOM\n')
        for num, key in enumerate(names.items()):
            key, item = key
            f1.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(num+1, item, str(positions[key][0]),
                                                   str(positions[key][1]),
                                                   str(positions[key][2]), item,
                                                                  # filter(lambda x: x.isalpha(), item),
                                                                       '1', 'LIG1', '0'))
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
###########################################mol2_files##############
    ln = coordinates_to_mm_basis_notation(xyz_names_bonds('Caffein.mol2'), valid_length=0.5, save_length=True)
    # ln = coordinates_to_notation(xyz_names('benzene.mol2'), save_length=True)
    print("coords_in_notation", ln)
    # print(ln[2],'2')
    print(to_one_way_bonds(ln[0]))
    print(ln)
    d31, nms = dimensional_length_unique_basis(ln)

    print(to_one_way_bonds(ln[0]))
    write_mol_file('My_caffein.mol2', nms, d31, bonds=to_one_way_bonds(ln[0]), attrs=ln[2])