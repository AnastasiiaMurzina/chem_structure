#!bin/bash
import sys, os, subprocess, shutil, tempfile
import numpy as np
from itertools import product
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file

if __name__ == '__main__':
    atom_name = "Caffein"
    '''
    Methods: 'first' - perpendicular rotations;
    'incline' - non-perpendicular (also possible choose second axis)
    'ten' - fixed ten rotation position (radial rotation)'''
    with open(atom_name+'compare_methods.txt', 'w') as comp:
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
            subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file1], stdout=subprocess.DEVNULL)
            with open(sim_file1, 'r') as f:
                f.readline()
                coeffs = f.readline().split()
                coeffs = coeffs[2:6] + coeffs[8:]
            comp.write('\n n_y = ' + str(i) + ' n_z = ' + str(j) + '\n')
            comp.write('\t'.join(coeffs))
        comp.write('\n incline method: \n')
        mth = 'incline'
        for fr, sr in product([1, 2], range(6, 11)):
            comp.write('\n incline '+str(fr)+' '+str(sr)+'\n')
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
                subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file], stdout=subprocess.DEVNULL)
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
        subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file], stdout=subprocess.DEVNULL)
        with open(sim_file, 'r') as f:
            f.readline()
            coeffs = f.readline().split()
            coeffs = coeffs[2:6] + coeffs[8:]
            comp.write('\t'.join(coeffs))
            # print(coeffs)
            # coeffs = np.array(list(map(float, coeffs)))
            # print(np.linalg.norm(coeffs))
        shutil.rmtree(tmpdir)
