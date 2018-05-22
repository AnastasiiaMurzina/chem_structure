#!bin/bash
import os, subprocess, shutil, tempfile
import numpy as np
from penta_with_rotate import get_penta_points, find_section,\
     Ry, Rz, d_ver_angle, d_hor_angle, step_rot
from itertools import product
from mol2_chain import bonds_to_one_way_dict, prepare_bonds, to_two_ways_bond2
import copy
from mpl_toolkits.mplot3d import axes3d
from mol2_chain import atoms_and_bonds, xyz_names_bonds, bonds_of_paired, write_mol2_file
import matplotlib.pyplot as plt
from matplotlib import cm
pp = get_penta_points([0,0,0])


def rotate_ten_vars(point, i1, left=-1, right=1):
    b = list(product([left, 0, right], repeat=2))
    operator = Rz(d_hor_angle*b[i1][0]).dot(Ry(d_ver_angle*b[i1][1]))
    pp_ = [i.dot(operator) for i in pp]
    for i in pp:
        pp_.append(i)
    if isinstance(point, (np.ndarray)):
        return point.dot(operator)
    return [i.dot(operator) for i in point]


def find_section(p0, p1, basis0=np.zeros(2), let_accurance=step_rot, all_posibility=False, left=-1,right=1):
    '''
    :param p0: this point has already basis
    :param p1: not important basis of p1
    :param basis0: basis of p0-center atom
    :return: section of p0 atom in which there's p1
    '''
    vec = p1 - p0
    vec = np.array([i/np.linalg.norm(vec) for i in vec])
    pp_ = rotate_ten_vars(pp, basis0, left=left, right=right)
    if all_posibility:
        return min([[np.linalg.norm(ppx - vec), ix] for ix, ppx in enumerate(pp_)])[1]
    while True:
        for num, i in enumerate(pp_):
            if np.linalg.norm(i - vec) <= let_accurance:
                return num
        let_accurance *= 1.1


def find_basis(point, connected, left=-1, right=1):
    '''
    :param point: point for search basis
    :param connected: atoms which have bonds with point
    :return: basis for point (y, z) by min of max different between point and center of section
    '''
    diffs = []
    for j in range(9):
        diff = []
        for i in connected:
            v = i - point
            v /= np.linalg.norm(v)
            diff.append(min([np.linalg.norm(v - ppx) for ppx in
                             rotate_ten_vars(pp, j, left=left, right=right)]))
        diffs.append([max(diff), j])
        # return min(diff)[1]
    return min(diffs)[1]


def mol2_to_notation(info_from_file, left=-1, right=1):
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
        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position(), basis0=basis, left=left, right=right)]
                                     for i in connected]), basis]})
        atoms[key].set_orientation(basis)

        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position(), basis0=basis, left=left, right=right)]
                                     for i in connected]), basis]})
    for key, item in bonds.items():
        for i in range(len(item)):
            bonds[key][i].insert(1, np.linalg.norm(atoms[key].position()-atoms[item[i][0]].position()))
    return notation, bonds


def dimensional_structure(notation, left=-1, right=1):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    with length
    :return: xyz-info
    '''
    dim_structure = {1: np.array([0, 0, 0])}#, np.array([0, 0])]}
    bonds_l, lengths = notation
    p = bonds_l[1]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, 1)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        for i in bonds:  # build bonds for cur_key atom
            if not (i[0] in dim_structure):  # if we don't have position:
                coord = rotate_ten_vars(pp[i[1]], basis, left=left, right=right)*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                                                 else lengths[(i[0], cur_key)][0]) \
                                + dim_structure[cur_key]

                dim_structure.update({i[0]: coord})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure


def symmetrical_errors(atom_name):
    '''
    :param atom_name: begin of mol2-file name for reading structure
    :return: graphic depend correctly coefficient of parameter ten-rotated basis method
    '''
    x = np.linspace(-1, 0, 50)
    y = []
    tmpdir = tempfile.mkdtemp()

    d = os.getcwd()
    compared_molecule = os.path.join(d, atom_name + '.mol2')
    atoms_info = atoms_and_bonds(compared_molecule)
    for left in x:
        sim_file = os.path.join(tmpdir, "sim_ten" + str(left) + ".txt")
        ln = mol2_to_notation(xyz_names_bonds(compared_molecule), left=left, right=-left)
        paired = bonds_of_paired(ln[1])
        dim_structure = dimensional_structure([ln[0], paired], left=left, right=-left)
        compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(atom_name))
        write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
        subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                        stdout=subprocess.DEVNULL)
        with open(sim_file, 'r') as f:
            f.readline()
            coeffs = f.readline().split()
            coeffs = coeffs[2:6] + coeffs[8:]
            # print(coeffs)
            coeffs = np.array(list(map(float, coeffs)))
            y.append(np.linalg.norm(coeffs))
            # print(np.linalg.norm(coeffs))
    shutil.rmtree(tmpdir)
    print(max(y))
    plt.plot(x, y)
    plt.show()


def asymmetrical_errors(atom_name, save=False):
    '''
    :param atom_name: begin of mol2-file name for reading structure
    :param save: save in pdf - boolean
    :return: graphic depend correctly coefficient of parameter ten-rotated basis method
    '''
    x = np.linspace(-1, 0, 5)
    y = np.linspace(0, 1, 5)
    z = []
    tmpdir = tempfile.mkdtemp()
    d = os.getcwd()
    compared_molecule = os.path.join(d, atom_name + '.mol2')
    atoms_info = atoms_and_bonds(compared_molecule)
    for left in x:
        z1 = []
        for right in y:
            sim_file = os.path.join(tmpdir, "sim_ten" + str(left) + str(right) + ".txt")
            ln = mol2_to_notation(xyz_names_bonds(compared_molecule), left=left, right=right)
            paired = bonds_of_paired(ln[1])
            dim_structure = dimensional_structure([ln[0], paired], left=left, right=right)
            compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(atom_name))
            write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
            subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                            stdout=subprocess.DEVNULL)
            with open(sim_file, 'r') as f:
                f.readline()
                coeffs = f.readline().split()
                coeffs = coeffs[2:6] + coeffs[8:]
                coeffs = np.array(list(map(float, coeffs)))
                z1.append(np.linalg.norm(coeffs))
        z.append(z1)
    shutil.rmtree(tmpdir)
    z = np.array(z)
    x, y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                    linewidth=0, antialiased=False)
    if save:
        plt.savefig(atom_name, format='pdf')
    else:
        plt.show()


if __name__ == '__main__':
    atom_name = "solution_neural_Ru_eaxial"
    '''
    Methods: 'first' - perpendicular rotations;
    'incline' - non-perpendicular (also possible choose second axis)
    'ten' - fixed ten rotation position (radial rotation)'''
    # mth = 'ten'
    asymmetrical_errors(atom_name, save=False)
