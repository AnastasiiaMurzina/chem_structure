import copy
import os
from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np

from mol_api import Molecule, compare_structers
from quadro_with_rotate_class import Spherical_divider
from searcher_in_space import Equation_system, apply_solution, length_xyz_interaction, genetic_to_the_aim, read_report

def reject_chjecker(reject_count):
    return reject_count > 20 and np.random.random() < 1-2**(-reject_count//10)

def approx_genetic_to_the_aim(reactant: Molecule, product: Molecule, system, solution,
                             file_log='reaction_report', prob_of_approx=0.99, trasher='trasher',
                              alpha=10):
    """
    :return: energies_path and length of sections changes
    """
    exact_line = '##\n'
    approx_line = '#?\n'

    react = copy.deepcopy(reactant)
    counters = SimpleNamespace(exist_eq_counter=0, can_solve_counter=0, approx_counter=0, smop_counter=0)

    def apply_change():
        mutant.refresh_dimensional()
        names = [i.name for _, i in mutant.atoms.items()]
        interaction = length_xyz_interaction((names, mutant.to_positions_array()))
        appr_flag = False
        exist_eq = system.find_equation_in_system(interaction, get_solution=True)
        if exist_eq[0] != -1:
            en = exist_eq[1]
            counters.exist_eq_counter += 1
        elif system.can_i_solve_it(interaction):
            counters.can_solve_counter += 1
            if np.random.random() < prob_of_approx:
                counters.approx_counter += 1
                en = apply_solution(solution, interaction)
                appr_flag = True
        if not appr_flag:
            en = mutant.get_energy()
            counters.smop_counter += 1
            mutant.to_xyz(trasher, title=str(en), mode='a')
            with open(trasher, 'a') as f_e: f_e.write(exact_line)
        if en is None:
            return False
        delta = en - en_path[-1]
        flag = False
        if np.random.random() < np.exp(-delta/alpha):
            flag = True
            en_path.append(en)
            mutant.to_xyz(file_log, title=str(en), mode='a')
            with open(file_log, 'a') as f_w1:
                f_w1.write(approx_line if appr_flag else exact_line)
        return flag

    d = react.notation.diff(product.notation)
    en_path = [react.get_energy()]
    react.to_xyz(trasher, title=str(en_path[-1]), mode='a')
    with open(trasher, 'a') as f_e: f_e.write(exact_line)
    reject_counter, all_rejects = 0, 0
    while d[2] > 0.1 and d[1] != 0:
        mutant = copy.deepcopy(react)
        if np.random.random() < 0.5:
            mutant.notation.s_change_step(product.notation)
        else:
            mutant.notation.l_change_step(product.notation)
        if apply_change():
            react = copy.deepcopy(mutant)
            reject_counter = 0
        else:
            reject_counter += 1
            all_rejects += 1
            if reject_chjecker(reject_counter):
                return (-1, counters.exist_eq_counter, counters.can_solve_counter, counters.approx_counter, counters.smop_counter, all_rejects)
        d = react.notation.diff(product.notation)
    while mutant.notation.l_change_step(product.notation) != -1 or mutant.notation.s_change_step(product.notation) != -1:
        if apply_change():
            react = copy.deepcopy(mutant)
            reject_counter = 0
        else:
            reject_counter += 1
            all_rejects += 1
            if reject_chjecker(reject_counter):
                return (-1, counters.exist_eq_counter, counters.can_solve_counter, counters.approx_counter, counters.smop_counter, all_rejects)
        mutant = copy.deepcopy(react)
    return (en_path, counters.exist_eq_counter, counters.can_solve_counter, counters.approx_counter, counters.smop_counter, all_rejects)


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

    learn = 'all_eqs_1'
    test = 'all_eqs_2'
    n = 13
    system = Equation_system()
    import pickle, random
    cache_file = "divider%d.pkl" % n
    if not os.path.exists(cache_file):
        divider = Spherical_divider(n)
        with open(cache_file, 'wb') as outp:
            pickle.dump(divider, outp, pickle.HIGHEST_PROTOCOL)
    else:
        with open(cache_file, 'rb') as inp:
            divider = pickle.load(inp)
    system.from_reports([learn])
    ss = system.solve()
    pred = []
    ens, poss = read_report(test, xyz=True)
    to_names = poss.pop(0)
    names = [i[0] for i in to_names]
    atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
    interaction = length_xyz_interaction((names, atoms))
    pred.append(apply_solution(ss, interaction))
    ix = 0
    popper = []
    print(len(ens))
    for item, energy in zip(poss, ens):
        atoms = [np.array(list(map(float, i[1::]))) for i in item]
        interaction = length_xyz_interaction((names, atoms))
        if isinstance(interaction, dict) and system.can_i_solve_it(interaction):
            pred.append(apply_solution(ss, interaction))
            if abs(pred[-1]-ens[ix]) > 2:
                pred.pop(-1)
                popper.append(ix)
        else:
            popper.append(ix)
        ix +=1

    for i in popper[::-1]:
        ens.pop(i)
    print(len(system.variables), len(system.energy))
    print(np.corrcoef(np.array(ens), np.array(pred)))
    print(len(ens))
    print(np.mean(np.array([abs(i) for i in np.array(ens)-np.array(pred)])))
    plt.scatter(ens, pred)
    plt.xlabel('Mopac energy, eV')
    plt.ylabel('Approximate energy, eV')
    plt.show()

    # print('mopac calc2')
    # print("Exist equation\tCould solve approximately\t Solve approximately\t Mopac solved\tNumber of rejects")
    # while success_paths < 30:
    #     random.seed(all_paths+100+random_bias)
    #     tstep.append(-time())
    #     if len(system.energy) != 0:
    #         ss = system.solve()
    #     file_log = "{0}{1}".format(saver, str(success_paths))
    #     path, c_eq, c_can, c_calc, m_calc, rejs = approx_genetic_to_the_aim(ln, pr, system, ss, file_log=file_log, trasher=trasher, alpha=30.)
    #     all_paths += 1
    #     tstep[-1] += time()
    #     system = Equation_system()
    #     system.from_reports([trasher])
    #     if isinstance(path, int):
    #         try:
    #             os.remove(file_log)
    #         except FileNotFoundError:
    #             pass
    #         print('Path {}.Failed path'.format(str(all_paths)))
    #         print('{}\t{}\t{}\t{}\t{}'.format(str(c_eq), str(c_can), str(c_calc), str(m_calc), str(rejs)))
    #         continue
    #
    #     print('Successful {}th path'.format(str(success_paths)))
    #     print('{}\t{}\t{}\t{}\t{}'.format(str(c_eq), str(c_can), str(c_calc), str(m_calc), str(rejs)))
    #     ts = max(path)
    #     print('TS energy is ', str(ts))
    #     tss.append(ts)
    #     success_paths += 1
    #
    # print(tstep)
    # print(tss)




