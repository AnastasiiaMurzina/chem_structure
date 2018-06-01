import sys, os, subprocess, shutil, tempfile
import numpy as np
from itertools import product
import pandas as pd
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file
from test_with_all_structure_rotation import get_rotated

if __name__ == '__main__':
    molecules = 'Aniline', 'Caffein',  'Ethanol', 'Naphthalene', 'Water'#, 'solution_neural_Ru_eaxial'#, 'vacuum_cation_singlet_Fe_full', 'Water'
    # molecules = 'solution_neural_Ru_eaxial', #, 'vacuum_cation_singlet_Fe_full'
    errors = {}
    d = os.getcwd()
    tmpdir = tempfile.mkdtemp()
    sim_file1 = os.path.join(tmpdir, "sim_f.txt")

    for molecule in molecules:
        compared_molecule = os.path.join(d, molecule + '.mol2')
        compared_molecule2 = os.path.join(d, "My_" + molecule + '.mol2')
        rotated_name = 'Rot_{}.mol2'.format(molecule)
        compared_molecule_rotated = os.path.join(d, rotated_name)
        get_rotated(molecule)

        def comparator(compared_molecule=compared_molecule):
            atoms_info = atoms_and_bonds(compared_molecule)
            ln = mol2_to_notation(xyz_names_bonds(compared_molecule), n_y=i, n_z=j, method=mth)
            paired = bonds_of_paired(ln[1])
            dim_structure = dimensional_structure([ln[0], paired], n_y=i, n_z=j, method=mth)
            compared_molecule2 = os.path.join(d, 'My_{}.mol2'.format(molecule))
            write_mol2_file(compared_molecule2, atoms_info, dim_structure, bonds=paired)


        def shaep_comparator(compared_molecule2=compared_molecule2):
            subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule2, compared_molecule, sim_file1],
                            stdout=subprocess.DEVNULL)
            with open(sim_file1, 'r') as f:
                f.readline()
                coeffs = f.readline().split()
                coeffs = coeffs[2:6] + coeffs[8:]

        mth = 'first'
        # for i, j in product(range(3, 10), repeat=2):
        i, j = 4, 4
        error, coeffs = [], []
        comparator()
        comparator(rotated_name)
        shaep_comparator(compared_molecule2=compared_molecule_rotated)
        error.append(np.linalg.norm(np.array(coeffs)))
        errors.update({molecule: [min(error)]})
        # ('\n incline method: \n')
        print(mth, molecule)
        mth = 'incline'
        error, coeffs = [], []
        fr, sr = 1, 7
        # for fr, sr in product([1, 2], range(6, 11)):
        #     for i, j in product(range(3, 10), repeat=2):
        # i, j = 4, 4
        comparator()
        comparator(compared_molecule=compared_molecule_rotated)
        shaep_comparator(compared_molecule2=compared_molecule_rotated)
        error.append(np.linalg.norm(np.array(coeffs)))
        errors.update({molecule: [min(error)]})
        #('\n Ten-fixed\n')
        print(mth, molecule)
        mth = 'ten'
        i, j = None, None
        coeffs = []
        comparator()
        comparator(rotated_name)
        shaep_comparator(compared_molecule2=compared_molecule_rotated)
        error.append(np.linalg.norm(np.array(coeffs)))
        errors.update({molecule: [min(error)]})
        print(mth, molecule)

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