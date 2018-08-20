import sys, os, subprocess, shutil, tempfile
from itertools import combinations, permutations, combinations_with_replacement
from copy import deepcopy
import numpy as np
from random import random
from many_mols import xyz_names_bonds, mol2_to_notation, molecular_divider, get_notation_many_mols, places_for_zero_bonds, to_two_ways_bond, unpack_with_zero_bonds, zero_connecter, write_mol2_file
from mopac_worker import writeInputFile, get_energy, mopacOut_to_xyz

FNULL = open(os.devnull, 'w')

if __name__ == '__main__':
    mols_dir = './intermediates'
    # subprocess.call(['ls', mols_dir])
    options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
               'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7', 'TS': 'TRANS'}
    tmpdir = tempfile.mkdtemp()
    name = 'INT_3a'
    # ini_mol_xyz = os.path.join(mols_dir, name+'.xyz')
    # ini_mol_mol2 = os.path.join(tmpdir, name+'.mol2')
    ini_mol_mol2 = './intermediates/3a.mol2'
    # ini_mol_mol2 = './TS/2b-2c/TS-2b-2c_step_1.mol2'
    method_compress = 'ten'
    # subprocess.call(['babel', '-ixyz', ini_mol_xyz, '-omol2', ini_mol_mol2])

    # subprocess.call(['babel', ])

    bs, ass = xyz_names_bonds(ini_mol_mol2)
    div_atoms, ligs_as, lig_bonds = molecular_divider(ass, bs)
    atoms_notation, bonds_notation = get_notation_many_mols(ass, bs, method=method_compress)
    needs_to_zero_discribe = places_for_zero_bonds(ass, bs, method=method_compress)

    zero_bonds = zero_connecter(ass, needs_to_zero_discribe, method=method_compress)
    # print('before all: ', zero_bonds)
    zero_bonds_1 = deepcopy(zero_bonds)
    # print(zero_bonds_1)
    # for counter in range(1000000):
    for counter in range(1):
        # for k, i in zero_bonds_1.items():
        #     for ix, j in enumerate(i):
        #         zero_bonds_1[k][ix][0][0] *= 1.25
        #         zero_bonds_1[k][ix][-1] = int(random()*10)

    # for i in permutations(range(10), len(zero_bonds)):
    #     for k, j in zero_bonds_1.items():
    #         for ix, sec in enumerate(zero_bonds_1[k]):
    #             zero_bonds_1[k][ix][-1] = i[k-1]
    #     print(zero_bonds)

        code_name = str(counter)
        coords_1 = unpack_with_zero_bonds(deepcopy(atoms_notation), deepcopy(bonds_notation), deepcopy(zero_bonds_1), method=method_compress)
        looking_for_ts_mol2 = os.path.join(tmpdir, name+'_random_search_'+code_name+'.mol2')
        write_mol2_file(looking_for_ts_mol2, ass, coords_1, to_two_ways_bond(bs, with_attr=True))
        looking_for_ts_xyz = os.path.join(tmpdir, name+'_random_search_'+code_name+'.xyz')
        looking_for_ts_mop = os.path.join(tmpdir, name+'_random_search_'+code_name+'.mop')
        subprocess.call(['babel', '-imol2', looking_for_ts_mol2, '-oxyz', looking_for_ts_xyz])
        writeInputFile(options, coords_1, ass, looking_for_ts_mop)
        mopac_opt_file_name = os.path.join(tmpdir, name+'opt_random_search_'+code_name+'.xyz')
        mopac_out_file = os.path.join(tmpdir, name + '_random_search_' + code_name)
        subprocess.call(['/opt/mopac/run_script.sh', looking_for_ts_mop], stdout=FNULL)
        # print(get_energy(mopac_out_file+'.out'))
        # mopacOut_to_xyz(mopac_out_file, mopac_opt_file_name)

    coords = unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds, method=method_compress)

    # looking_for_ts_mol2 = os.path.join(tmpdir, name+'_random_search.mol2')
    # write_mol2_file(looking_for_ts_mol2, ass, coords, to_two_ways_bond(bs, with_attr=True))
    ##################################################

    # shutil.rmtree(tmpdir)


