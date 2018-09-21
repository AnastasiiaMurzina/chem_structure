#!bin/bash
import sys, os, subprocess, shutil, tempfile
from icosahedron_with_rotate import *
from penta_with_rotate import *
from itertools import product
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import numpy.random
from numpy import arccos, arctan, arctan2

color_map = dict(enumerate(["r", "g", "b", "peachpuff", "fuchsia", 'c', 'm', 'y', 'k', 'burlywood', 'chartreuse',
                                'olive', 'maroon', 'peru', 'sienna', 'seagreen', 'hotpink', 'lime', 'salmon', 'gray']))

def cartesian_to_spherical(vector):
    r_xy = vector[0] ** 2 + vector[1] ** 2
    theta = arctan2(vector[1], vector[0])
    phi = arctan2(vector[2], r_xy ** 0.5)
    return theta, phi

def smth_before_script():
    atom_name = "Caffein"
    '''
    Methods: 'first' - perpendicular rotations;
    'incline' - non-perpendicular (also possible choose second axis)
    'ten' - fixed ten rotation position (radial rotation)'''
    with open(atom_name + 'compare_methods.txt', 'w') as comp:
        mth = 'first'
        comp.write('first method:\n')
        tmpdir = tempfile.mkdtemp()
        for i, j in product(range(1, 10), repeat=2):
            sim_file1 = os.path.join(tmpdir, "sim_f_" + str(i) + str(j) + '.txt')
            d = os.getcwd()
            compared_molecule = os.path.join(d, atom_name + '.mol2')
            atoms_info = atoms_and_bonds(compared_molecule)
            ln = mol2_to_notation(xyz_names_bonds(compared_molecule), n_y=i, n_z=j, method=mth)
            paired = bonds_of_paired(ln[1])
            dim_structure = dimensional_structure([ln[0], paired], n_y=i, n_z=j, method=mth)
            compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(atom_name))
            write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
            subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file1],
                            stdout=subprocess.DEVNULL)
            with open(sim_file1, 'r') as f:
                f.readline()
                coeffs = f.readline().split()
                coeffs = coeffs[2:6] + coeffs[8:]
            comp.write('\n n_y = ' + str(i) + ' n_z = ' + str(j) + '\n')
            comp.write('\t'.join(coeffs))
        comp.write('\n incline method: \n')
        mth = 'incline'
        for fr, sr in product([1, 2], range(6, 11)):
            comp.write('\n incline ' + str(fr) + ' ' + str(sr) + '\n')
            for i, j in product(range(1, 10), repeat=2):
                tmpdir = tempfile.mkdtemp()
                sim_file = os.path.join(tmpdir, "sim_" + str(i) + str(j) + '.txt')
                d = os.getcwd()
                compared_molecule = os.path.join(d, atom_name + '.mol2')
                atoms_info = atoms_and_bonds(compared_molecule)
                ln = mol2_to_notation(xyz_names_bonds(compared_molecule), n_y=i, n_z=j, method=mth)
                paired = bonds_of_paired(ln[1])
                dim_structure = dimensional_structure([ln[0], paired], n_y=i, n_z=j, method=mth)
                compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(atom_name))
                write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
                subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                                stdout=subprocess.DEVNULL)
                with open(sim_file, 'r') as f:
                    f.readline()
                    coeffs = f.readline().split()
                    coeffs = coeffs[2:6] + coeffs[8:]
                comp.write('\n n_y = ' + str(i) + ' n_z = ' + str(j) + '\n')
                comp.write('\t'.join(coeffs))
                # comparator.update({(i, j): list(map(float, coeffs))})
                # print(comparator)
                # m = 1
        #     for i, j in comparator.items():
        #         m = max(m, np.linalg.norm(j))
        #     print(m)
        comp.write('\n Ten-fixed\n')
        mth = 'ten'
        tmpdir = tempfile.mkdtemp()
        sim_file = os.path.join(tmpdir, "sim_ten.txt")
        d = os.getcwd()
        compared_molecule = os.path.join(d, atom_name + '.mol2')
        atoms_info = atoms_and_bonds(compared_molecule)
        ln = mol2_to_notation(xyz_names_bonds(compared_molecule), method=mth)
        paired = bonds_of_paired(ln[1])
        dim_structure = dimensional_structure([ln[0], paired], method=mth)
        compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(atom_name))
        write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
        subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                        stdout=subprocess.DEVNULL)
        with open(sim_file, 'r') as f:
            f.readline()
            coeffs = f.readline().split()
            coeffs = coeffs[2:6] + coeffs[8:]
            comp.write('\t'.join(coeffs))
            # print(coeffs)
            # coeffs = np.array(list(map(float, coeffs)))
            # print(np.linalg.norm(coeffs))
        shutil.rmtree(tmpdir)

def graph_ico_on_2d():
    for i in icos:
        theta, phi = cartesian_to_spherical(i)
        plt.scatter(theta, phi)
    plt.xlim(-pi, pi)
    plt.ylim(-pi/2, pi/2)
    plt.show()

def graphic_first_ico():
    for j in product(range(9), repeat=2):
        icos_ = rotate_by_basis(icos, j[0], j[1])
        for num, i in enumerate(icos_):
            theta, phi = cartesian_to_spherical(i)
            plt.scatter(theta, phi)
    plt.show()

def graphic_incline_ico():
    for j in product(range(9), repeat=2):
        icos_ = rotate_non_perpendicular(icos, j[0], j[1])
        for i in icos_:
            theta = arccos(i[2])
            phi = arctan(i[1] / i[0])
            plt.scatter(theta, phi)
    plt.show()

def graphic_ten_ico():
    for j in range(9):
        icos_ = rotate_ten_vars(icos, j)
        for i in icos_:
            theta = arccos(i[2])
            phi = arctan(i[1] / i[0])
            plt.scatter(theta, phi)
    plt.show()

def graphic_first_ico_color_by_num(n_y=9, n_z=9, **kwargs):
    for j in product(range(n_y), range(n_z)):
        icos_ = rotate_by_basis(icos, j[0], j[1], n_y=n_y, n_z=n_z, kwargs=kwargs)
        for num, i in enumerate(icos_):
            theta, phi = cartesian_to_spherical(i)
            plt.scatter(theta, phi, c=color_map[num%20])
    plt.show()

def graphic_incline_ico_color_by_num(n_y=9, n_z=9,**kwargs):
    for j in product(range(n_y), range(n_z)):
        icos_ = rotate_non_perpendicular(icos, j[0], j[1], n_y=n_y, n_z=n_z, **kwargs)
        for num, i in enumerate(icos_):
            theta, phi = cartesian_to_spherical(i)
            plt.scatter(theta, phi, c=color_map[num%20])
    plt.show()

def graphic_ten_ico_color_by_num(r=0.6):
    for j in range(9):
        icos_ = rotate_ten_vars(icos, j, r=r)
        for num, i in enumerate(icos_):
            theta, phi = cartesian_to_spherical(i)
            plt.scatter(theta, phi, c=color_map[num%20])
    plt.show()

def graphic_vars_ico_color_by_num(r=0.6, ns=3):
    for j in range(9):
        icos_ = rotate_radius_vars(icos, j, r=r, ns=ns)
        for num, i in enumerate(icos_):
            theta, phi = cartesian_to_spherical(i)
            plt.scatter(theta, phi, c=color_map[num%20])
    plt.show()

def graphic_errors(mth='first', **kwargs):
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    n = 8000
    method = mth
    n_y, n_z = kwargs.get('n_y', 4), kwargs.get('n_z', 4)
    fr, sr = kwargs.get('fr', 1), kwargs.get('sr', 10)
    color_map = plt.get_cmap('cool')
    average_error = 0
    fig = plt.figure()
    ax1 = plt.subplot(111)
    mx_error = 0
    for _ in range(n):
        x, y, z = numpy.random.normal(size=3)
        r = (x**2+y**2+z**2)**0.5
        x /= r
        y /= r
        z /= r
        theta, phi = cartesian_to_spherical(np.array([x,y,z]))
        cur_err = get_error_of_point(np.array([x, y, z]), method=method, n_y=n_y, n_z=n_z)
        d_err = mx_error - average_error
        sc = ax1.scatter(theta, phi, c=cur_err, cmap=color_map, vmin=average_error-d_err/2, vmax=average_error+d_err/2)
        average_error += cur_err/n
        mx_error = max(mx_error, cur_err)
    ax1.set_title(method+', icosahedron, n_y='+str(n_y)+', n_z='+str(n_z))
    cb = fig.colorbar(sc, orientation='horizontal', drawedges=False)
    cb.set_label('Discrete errors, average='+str(format(average_error, '.5g'))+', Maximal = ' + str(format(mx_error, '.5g')), fontsize=14)
    plt.show()

if __name__ == '__main__':
    
    # graphic_errors(mth='first', n_y=2, n_z=2)
    # graphic_errors(mth='first', n_y=4, n_z=4)
    # graphic_errors(mth='first', n_y=10, n_z=10)
    # graphic_errors(mth='first', n_y=20, n_z=20)

    # graphic_errors(mth='incline', n_y=2, n_z=2, fr=1, sr=7)
    # graphic_errors(mth='incline', n_y=2, n_z=2, fr=1, sr=8)
    # graphic_errors(mth='incline', n_y=2, n_z=2, fr=1, sr=9)
    # graphic_errors(mth='incline', n_y=2, n_z=2, fr=1, sr=10)

    # graphic_errors(mth='ten', r=0.75)


    # graph_ico_on_2d()

    # graphic_first_ico_color_by_num(n_y=9, n_z=9)
    # graphic_first_ico_color_by_num(n_y=4, n_z=7)
    # graphic_first_ico_color_by_num(n_y=4, n_z=4)
    # graphic_first_ico_color_by_num(n_y=20, n_z=20)
    # graphic_incline_ico_color_by_num(n_y=5, n_z=5)
    # graphic_incline_ico_color_by_num(n_y=20, n_z=20, fr=1,sr=7)
    # graphic_incline_ico_color_by_num(n_y=20, n_z=20, fr=1,sr=8)
    # graphic_incline_ico_color_by_num(n_y=20, n_z=20, fr=1,sr=9)
    # graphic_incline_ico_color_by_num(n_y=20, n_z=20, fr=1,sr=10)

    graphic_ten_ico_color_by_num(0.666667)
    # graphic_ten_ico_color_by_num(0.75)
    # graphic_ten_ico_color_by_num(0.9)
    # graphic_ten_ico_color_by_num(0.8)
    # graphic_ten_ico_color_by_num(0.65)
    # graphic_vars_ico_color_by_num(0.75, 10)
    # graphic_vars_ico_color_by_num(0.75, 200)



