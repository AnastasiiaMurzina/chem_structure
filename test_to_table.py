import sys, os, subprocess, shutil, tempfile
import numpy as np
from itertools import product
import pandas as pd
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file


if __name__ == '__main__':
    # molecules = 'Aniline', 'Caffein',  'Ethanol', 'Naphthalene', 'Water'#, 'solution_neural_Ru_eaxial'#, 'vacuum_cation_singlet_Fe_full', 'Water'
    molecules = 'solution_neural_Ru_eaxial', #, 'vacuum_cation_singlet_Fe_full'
    errors = {}
    for molecule in molecules:
        mth = 'first'
        tmpdir = tempfile.mkdtemp()
        error = []
        for i, j in product(range(1, 10), repeat=2):
            sim_file1 = os.path.join(tmpdir, "sim_f_" + str(i) + str(j) + '.txt')
            d = os.getcwd()
            compared_molecule = os.path.join(d, molecule + '.mol2')
            atoms_info = atoms_and_bonds(compared_molecule)
            ln = mol2_to_notation(xyz_names_bonds(compared_molecule), n_y=i, n_z=j, method=mth)
            paired = bonds_of_paired(ln[1])
            dim_structure = dimensional_structure([ln[0], paired], n_y=i, n_z=j, method=mth)
            compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(molecule))
            write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
            subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file1],
                            stdout=subprocess.DEVNULL)
            with open(sim_file1, 'r') as f:
                f.readline()
                coeffs = f.readline().split()
                coeffs = coeffs[2:6] + coeffs[8:]
            error.append(np.linalg.norm(np.array(coeffs)))
        errors.update({molecule: [min(error)]})
        # ('\n incline method: \n')
        mth = 'incline'
        error = []
        for fr, sr in product([1, 2], range(6, 11)):
            for i, j in product(range(1, 10), repeat=2):
                tmpdir = tempfile.mkdtemp()
                sim_file = os.path.join(tmpdir, "sim_" + str(i) + str(j) + '.txt')
                d = os.getcwd()
                compared_molecule = os.path.join(d, molecule + '.mol2')
                atoms_info = atoms_and_bonds(compared_molecule)
                ln = mol2_to_notation(xyz_names_bonds(compared_molecule), n_y=i, n_z=j, method=mth)
                paired = bonds_of_paired(ln[1])
                dim_structure = dimensional_structure([ln[0], paired], n_y=i, n_z=j, method=mth)
                compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(molecule))
                write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
                subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                                stdout=subprocess.DEVNULL)
                with open(sim_file, 'r') as f:
                    f.readline()
                    coeffs = f.readline().split()
                    coeffs = coeffs[2:6] + coeffs[8:]
                    error.append(np.linalg.norm(np.array(coeffs)))
        errors[molecule].append(min(error))
        #('\n Ten-fixed\n')
        mth = 'ten'
        error = []
        tmpdir = tempfile.mkdtemp()
        sim_file = os.path.join(tmpdir, "sim_ten.txt")
        d = os.getcwd()
        compared_molecule = os.path.join(d, molecule + '.mol2')
        atoms_info = atoms_and_bonds(compared_molecule)
        ln = mol2_to_notation(xyz_names_bonds(compared_molecule), method=mth)
        paired = bonds_of_paired(ln[1])
        dim_structure = dimensional_structure([ln[0], paired], method=mth)
        compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(molecule))
        write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)
        subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule, compared_molecule2, sim_file],
                        stdout=subprocess.DEVNULL)
        with open(sim_file, 'r') as f:
            f.readline()
            coeffs = f.readline().split()
            coeffs = coeffs[2:6] + coeffs[8:]
        errors[molecule].append(np.linalg.norm(np.array(coeffs)))
        shutil.rmtree(tmpdir)
    # molecule = molecules[-1]
    # bs, ass = xyz_names_bonds(molecule + '.mol2')
    # # div_atoms, ligs_as, lig_bonds = molecular_divider(ass, bs)
    # from many_mols import  get_notation_many_mols, places_for_zero_bonds, zero_connecter, unpack_with_zero_bonds, write_mol2_file
    # atoms_notation, bonds_notation = get_notation_many_mols(ass, bs)
    #
    # needs_to_zero_discribe = places_for_zero_bonds(ass, bs)
    #
    # zero_bonds = zero_connecter(ass, needs_to_zero_discribe)
    # coords = unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds)
    # from mol2_chain import to_two_ways_bond
    # write_mol2_file('long.mol2', ass, coords, to_two_ways_bond(bs, with_attr=True))
    # print(errors)

    df = pd.DataFrame(data=errors)
    df.index = ['first', 'incline', 'ten']
    print(df)