import copy
from tempfile import mkdtemp
from time import time

import matplotlib.pyplot as plt
import numpy as np

from layouter import dimensional_structure
from mol_api import Molecule, compare_structers
from quadro_with_rotate_class import Spherical_divider
from searcher_in_space import Equation_system, apply_solution, length_xyz_interaction, genetic_to_the_aim


def approx_genetic_to_the_aim(reactant: Molecule, product: Molecule, system, solution,
                             file_log='reaction_report', prob_of_approx=0.5, trasher='trasher'):
    """
    :return: energies_path and length of sections changes
    """
    exact_line = '##\n'
    approx_line = '#?\n'

    def apply_change():
        mutant.refresh_dimensional()
        names = [i.name for _, i in mutant.atoms.items()]
        interaction = length_xyz_interaction((names, mutant.to_positions_array()))
        appr_flag = False
        exist_eq = system.find_equation_in_system(interaction, get_solution=True)
        if exist_eq[0] != -1:
            en = exist_eq[1]
        elif system.can_i_solve_it(interaction) and np.random.random() < prob_of_approx:
            en = apply_solution(solution, interaction)
            appr_flag = True
        else:
            en = mutant.get_energy()
            mutant.to_xyz(trasher, title=str(en), mode='a')
            with open(trasher, 'a') as f_e: f_e.write(exact_line)
        if en is None:
            return False
        delta = en - path[-1]
        flag = False
        if np.random.random() < np.exp(-delta):
            flag = True
            path.append(en)
            mutant.to_xyz(file_log, title=str(en), mode='a')
            with open(file_log, 'a') as f_w1:
                f_w1.write(approx_line if appr_flag else exact_line)
        return flag

    d = reactant.notation.diff(product.notation)
    path = [reactant.get_energy()]
    reactant.to_xyz(trasher, title=str(path[-1]), mode='a')
    with open(trasher, 'a') as f_e: f_e.write(exact_line)
    while d[2] > 0.1 and d[1] != 0:
        mutant = copy.deepcopy(reactant)
        if np.random.random() < 0.5:
            mutant.notation.s_change_step(product.notation)
        else:
            mutant.notation.l_change_step(product.notation)
        if apply_change(): reactant = copy.deepcopy(mutant)
        d = reactant.notation.diff(product.notation)
    while mutant.notation.l_change_step(product.notation) != -1 or mutant.notation.s_change_step(product.notation) != -1:
        if apply_change(): reactant = copy.deepcopy(mutant)
        mutant = copy.deepcopy(reactant)
    return path


def read_approx_report(report_name):
    with open(report_name, 'r') as f:
        k = int(f.readline())
        ens = [float(f.readline())]
        accurance, poss = [], []
        flag = True
        while flag:
            mol = []
            for _ in range(k):
                mol.append(f.readline().split())
            accurance.append(not ('?' in f.readline()))
            flag = f.readline().isdigit()
            poss.append(mol)
            ens.append(float(f.readline()))
    return (ens, poss, accurance)


if __name__ == '__main__':

    def calc_to_the_aim_path(n):
        divider = Spherical_divider(n=n)
        reaction = 'mopac_example' # '3a->4' #
        # reaction = '3a->4' #
        # reaction = 'vanadii'
        if reaction == '3a->4':
            ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n, divider=divider)
            pr = Molecule('./prepared_mols2/4_opted.mol2', n=n, divider=divider)

        # elif reaction == 'vanadii':
        #     ln = Molecule('./vanadii/3a_singlet_opted.mol2', n=n)
        #     pr = Molecule('./vanadii/ts_3a_4a_opted.mol2', n=n)
        else:
            ln = Molecule('./ordered_mol2/js_exapmle_init.mol2', n=n)
            pr = Molecule('./ordered_mol2/js_exapmle_finish.mol2', n=n)

        # kk = 115
        start_energy = ln.get_energy()
        finish_energy = pr.get_energy()
        pr.refresh_dimensional()
        ms = genetic_to_the_aim(ln, pr, file_log=reaction + '_report_2')#+str(kk))
        ms.insert(0, start_energy)
        ms.append(finish_energy)
        print(max(ms))
        # print(kk, max(ms))

    # system = Equation_system()
    # # system.from_reports(['equations_mn_3a4_0', 'equations_mn_3a4_1'])
    # system.from_reports(['equations_mn_3a4_0', 'equations_mn_3a4_1'])
    # ss = system.solve()
    # print(len(ss))
    # print(len(system.variables))

    file_name = './ordered_mol2/js_exapmle_init.mol2'
    # file_name = './prepared_mols2/3a_opted.mol2'
    to_file = './ordered_mol2/js_exapmle_finish.mol2'
    # to_file = './prepared_mols2/4_opted.mol2'
    saver = 'approx_report_path0_'
    # saver = 'equations_mnw_3a4_4'
    trasher = 'all_eqs_0'
    n = 20
    system = Equation_system()
    divider = Spherical_divider(n)
    ss = {}
    tstep = []
    for ix in range(5):
        tstep.append(-time())
        if len(system.energy) != 0:
            ss = system.solve()
        ln = Molecule(file_name, divider=divider)
        pr = Molecule(to_file, divider=divider)
        approx_genetic_to_the_aim(ln, pr, system, ss, file_log="{0}{1}".format(saver, str(ix)),
                                        trasher=trasher)
        tstep[-1] += time()
        system = Equation_system()
        system.from_reports([trasher])

    print(tstep)




