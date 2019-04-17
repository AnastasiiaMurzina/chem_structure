import copy
import numpy as np
from numpy import arctan2, pi
from mol2_chain_q import atoms_and_bonds,  bonds_of_paired, xyz_names_bonds, Notation
from mol2_worker import xyz_names_bonds
from layouter import check, relaxing, dimensional_structure

if __name__ == '__main__':

    name = './sep_experiment/3a_small_opt'
    n_param = 6

    atoms_info = atoms_and_bonds(name + '.mol2')
    ln_small = Notation(n=n_param, info_from_file=xyz_names_bonds(name + '.mol2'))
    # print(ln.notation)
    # print(ln.bonds)
    # for i in ln.atoms:
    #     print(i.)
    paired = bonds_of_paired(ln_small.bonds)
    dim_structure = dimensional_structure(ln_small, relax=True)

