import copy
from tempfile import mkdtemp

import matplotlib.pyplot as plt
import numpy as np

from layouter import dimensional_structure
from mol_api import Molecule, compare_structers
from quadro_with_rotate_class import Spherical_divider


def searcher(substrat, product, zero_bonds_only=False,
             length_change=True, length_path=10):
    '''
    :param substrat: type: Molecule,
    :param product: type: Molecule,
    :param zero_bonds_only: not ready
    :return: path of rmsd
    '''
    # c1 = substrat.compare_with(product.to_positions_array())
    dim_structure = dimensional_structure(substrat.notation, relax=True)
    c1 = product.compare_with(np.array([i for _, i in dim_structure.items()]))
    paths = [c1]
    st = copy.deepcopy(substrat)
    s_comp = substrat.notation.diff(product.notation)[1]
    while len(paths) < length_path:
    # while paths[-1] > 0.25:

        ch = st.get_child(zero=zero_bonds_only, change_length=length_change)
        ch_dim = dimensional_structure(ch.notation, relax=True)
        compared = product.compare_with(np.array([i for _, i in ch_dim.items()]))
        _, s_comp_c, _ = st.notation.diff(product.notation)
        if s_comp_c < s_comp:
        # if 1-np.random.rand() < c1/compared:
            paths.append(compared)
            print('structure is ', ch_dim)
            del st
            st = copy.deepcopy(ch)
    return paths

def real_random_path(reactant, product, n=10000, write=False, file_log='random_search_report'):
    '''
    :param reactant, product: Molecule
    :param n: number of loops
    :return: show plot of rmsds
    '''
    def apply_change():
        mut.refresh_dimensional()
        # path.append(reactant.notation.get_heat())
        msds.append(compare_structers(mut.to_positions_array(), product.to_positions_array()))
        if write:
            mut.to_xyz(file_log, mode='a')
            with open(file_log, 'a') as f_w:
                f_w.write('################\n')
    with open(file_log, 'w') as f_w:
        f_w.write('\n\n\n\n')
    mut = copy.deepcopy(reactant)
    d = reactant.notation.diff(product.notation)
    msd = compare_structers(mut.to_positions_array(), product.to_positions_array())
    mut_pr = reactant.mutation()
    mut_pr.refresh_dimensional()
    d_pr = mut_pr.notation.diff(product.notation)
    msd_pr = compare_structers(mut_pr.to_positions_array(), product.to_positions_array())
    msds = [msd]
    apply_change()
    print('initial msd', msd, 'd', d)
    print('rmsd accept probability')
    for _ in range(n):
        # if d_pr[1] < d[1] or d_pr[2] < d[2] or np.random.random() < np.exp(-(d_pr[2]/d[2])):
        if d_pr[1] < d[1] or d_pr[2] < d[2] or msd_pr < msds[0]*1.5:
            mut, d, msd = copy.deepcopy(mut_pr), d_pr, msd_pr
            msds.append(msd)
            apply_change()
        mut_pr = mut.mutation()
        mut_pr.refresh_dimensional()
        d_pr = mut_pr.notation.diff(product.notation)
        msd_pr = compare_structers(mut_pr.to_positions_array(), product.to_positions_array())
    print('final msd', msd, 'd', d)
    if write:
        f = open(file_log, 'r+')
        f.seek(0, 0)
        f.write(str(len(msds)-1)+'\n')
        f.close()
    plt.plot(list(range(len(msds))), msds)
    plt.show()


def random_to_the_aim_search(reactant, product, write=False, file_log='reaction_report'): #TODO implement here add of bonds
    """
    :return: energies_path and length of sections changes
    """
    def apply_change():
        reactant.refresh_dimensional()
        # path.append(reactant.notation.get_heat())
        msds.append(compare_structers(reactant.to_positions_array(), product.to_positions_array()))
        if write:
            reactant.to_xyz(file_log, mode='a')
            with open(file_log, 'a') as f_w:
                f_w.write('################\n')
    d = reactant.notation.diff(product.notation)
    with open(file_log, 'w') as f_w:
        f_w.write(str(len(d))+'\n')
    # path = [reactant.notation.get_heat()]
    msds = []
    apply_change()
    while d[2] > 0.1 and d[1] != 0:
        if np.random.random() < 0.5:
            reactant.notation.s_change_step(product.notation)
        else:
            reactant.notation.l_change_step(product.notation)
        apply_change()
        d = reactant.notation.diff(product.notation)
    while reactant.notation.l_change_step(product.notation) != -1:
        apply_change()
    while reactant.notation.s_change_step(product.notation) != -1:
        apply_change()
    if False and write:
        f = open(file_log, 'r+')
        f.seek(0, 0)
        f.write(str(len(msds))+'\n')
        f.close()

    # return path, msds
    return msds


def genetic_to_the_aim(reactant: Molecule, product: Molecule, write=False,
                             file_log='reaction_report', trash_keeper=''):
    """
    :return: energies_path and length of sections changes
    """
    def apply_change():
        mutant.refresh_dimensional()
        en = mutant.get_energy()
        if en is None:
            return False
        delta = en - path[-1]
        flag = False
        if np.random.random() < np.exp(-delta):
            flag = True
            path.append(en)
            # print('difference accept', mutant.notation.l_diff(product.notation))
            if write:
                mutant.to_xyz(file_log, title=str(en), mode='a')
                with open(file_log, 'a') as f_w:
                    f_w.write('##\n')
        if trash_keeper != '':
            mutant.to_xyz(trash_keeper, title=str(en), mode='a')
            with open(trash_keeper, 'a') as f_w:
                f_w.write('##\n')
        return flag

    d = reactant.notation.diff(product.notation)
    path = [reactant.get_energy()]
    # msds = []
    # apply_change()
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
    return path #, msds

    #########two halves of random path#########
    # s = reactant.notation.s_change(product.notation, follow_energy=True, sh=True)
    # l = reactant.notation.l_change(product.notation, follow_energy=True, sh=True)
    # return s+l, len(s)
    ###########################################

def read_report(report_name, xyz=False):
    with open(report_name, 'r') as f:
        k = int(f.readline())
        ens = [float(f.readline())]
        if xyz: poss = []
        flag = True
        while flag:
            mol = []
            for _ in range(k):
                mol.append(f.readline().split())
            if xyz: poss.append(mol)
            f.readline()
            flag = f.readline().replace('\n','').isdigit()
            if flag: ens.append(float(f.readline()))
    return (ens, poss) if xyz else ens

def length_xyz_interaction(report_read, l_limit=0.5, g_limit=3.5, bias=True):
    names, atoms = report_read
    interaction = {tuple(['const', 'const', 0.0]): 1} if bias else {}
    for ix, i in enumerate(atoms):
        for jx, j in enumerate(atoms):
            if ix != jx:
                d = round(np.linalg.norm(i-j), 1)
                if l_limit > d:
                    return np.inf
                elif d < g_limit:
                    a1, a2 = names[ix], names[jx]
                    if a1 > a2: a1, a2 = a2, a1
                    cur_el = tuple([a1, a2, d])
                    if not interaction.get(cur_el):
                        interaction.update({cur_el: 0})
                    interaction[cur_el] += 1
    return interaction

class Equation_system:
    def __init__(self, first_linear=None, energy=None, bias=True):
        self.variables = set({})
        # self.variables = {tuple(['const', 'const', 0.0])} if bias else set([]) # probably if use this __init__
        self.equations = []
        self.energy = []
        if not (first_linear is None):
            self.equations.append(first_linear)
            self.variables.update(first_linear.keys())
        if not(energy is None):
            self.energy.append(energy)

    def to_matrix(self):
        x = np.zeros((len(self.energy), len(self.variables)))
        order_of_vars = []
        for ix_e, var in enumerate(self.variables):
            order_of_vars.append(var)
            for ix, i in enumerate(self.equations):
                x[ix][ix_e] = i.get(var, 0)
        return x, np.array(self.energy), order_of_vars

    def from_reports(self, files):
        for file in files:
            ens, poss = read_report(file, xyz=True)
            to_names = poss.pop(0)
            names = [i[0] for i in to_names]
            atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
            interaction = length_xyz_interaction((names, atoms))
            if isinstance(interaction, dict):
                self.push(interaction, ens.pop(0))
            for item, energy in zip(poss, ens):
                atoms = [np.array(list(map(float, i[1::]))) for i in item]
                interaction = length_xyz_interaction((names, atoms))
                if isinstance(interaction, dict):
                    self.push(interaction, energy)

    def push(self, linear, energy):
        self.equations.append(linear)
        self.variables.update(linear.keys())
        self.energy.append(energy)

    def check_system(self):
        return len(self.energy) > len(self.variables)

    def can_i_solve_it(self, variables):
        return len(variables.keys() - self.variables) == 0

    def solve(self):
        system, ens, order = self.to_matrix()
        try:
            return {i: j for i, j in zip(order, np.linalg.lstsq(system, ens, rcond=None)[0])}
        except np.linalg.LinAlgError:
            return -1

    def find_equation_in_system(self, interaction, get_solution=False):
        for i, e in enumerate(self.equations):
            if e == interaction:
                if get_solution:
                    return (i, self.energy[i])
                return i
        return (-1, -1) if get_solution else -1

def apply_solution(solution: dict, interacted: dict):
    en = 0
    for key, item in interacted.items():
        en += solution[key]*item
    return en


if __name__ == '__main__':
    def calc_to_the_aim_path(n):
        divider = Spherical_divider(n=n)
        reaction = '3a->4' #
        ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n, divider=divider)
        pr = Molecule('./prepared_mols2/4_opted.mol2', n=n, divider=divider)
        start_energy = ln.get_energy()
        finish_energy = pr.get_energy()
        pr.refresh_dimensional()
        ms = genetic_to_the_aim(ln, pr, write=True, file_log=reaction + '_report_2')#+str(kk))
        ms.insert(0, start_energy)
        ms.append(finish_energy)
        print(ms)
        print(max(ms))
        # print(kk, max(ms))
    print('15')
    calc_to_the_aim_path(15)

        # divider = Spherical_divider(n=n)
        # ln = Molecule(file_name, n=n, divider=divider)
        # # ln_init = copy.deepcopy(ln)
        # pr = Molecule(to_file, n=n, divider=divider)
        # lstl_solver = Equation_system()
        # lstl_solver.from_file(file_eqs)
        # ss = lstl_solver.solve()
        # d = ln.notation.diff(pr.notation)
        # # print(apply_solution(ss, ln.length_interaction()))
        # ln.refresh_dimensional()
        # mops = []
        # appr = []
        # while d[2] > 0.1 and d[1] != 0:
        #     if np.random.random() < 0.5:
        #         ln.notation.s_change_step(pr.notation)
        #     else:
        #         ln.notation.l_change_step(pr.notation)
        #     apply_change()
        #     d = ln.notation.diff(pr.notation)
        # while ln.notation.l_change_step(pr.notation) != -1:
        #     apply_change()
        # while ln.notation.s_change_step(pr.notation) != -1:
        #     apply_change()
        # print(mops)
        # print(appr)
        # print(np.corrcoef(np.array(mops), np.array(appr)))
        # plt.scatter(mops, appr)
        # plt.show()
        # print(file_name)


    def approx_report(report_file, ss, appr=[], mops=[]):

        ens, poss = read_report(report_file, xyz=True)
        to_names = poss.pop(0)
        names = [i[0] for i in to_names]
        atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
        interaction = length_xyz_interaction((names, atoms))
        if (not isinstance(interaction, float)) and system.can_i_solve_it(interaction):
            appr.append(apply_solution(ss, interaction))
            mops.append(ens.pop(0))

        for en, item in zip(ens, poss):
            atoms = [np.array(list(map(float, i[1::]))) for i in item]
            interaction = length_xyz_interaction((names, atoms))
            if (not isinstance(interaction, float)) and system.can_i_solve_it(interaction):
                appr.append(apply_solution(ss, interaction))
                mops.append(en)
        return mops, appr

    #
    # calc_to_the_aim_path(n)
    # system = Equation_system()
    # # system.from_reports(['equations_mn_3a4_0', 'equations_mn_3a4_1'])
    # system.from_reports(['equations_mn_3a4_0', 'equations_mn_3a4_1'])
    # ss = system.solve()
    # print(len(ss))
    # print(len(system.variables))

    # file_name = './ordered_mol2/js_exapmle_init.mol2'
    file_name = './prepared_mols2/3a_opted.mol2'
    # to_file = './ordered_mol2/js_exapmle_finish.mol2'
    to_file = './prepared_mols2/4_opted.mol2'
    saver = 'all_equations_mnw_3a4_rr'
    divider = Spherical_divider(n=10)

    # trasher = 'all_eqs_mnmn'
    n = 10
    system = Equation_system()
    system.from_reports(['all_equations_mnw_3a4_rr'])
    # system.from_reports(['equations_mopac_example_n0', 'equations_mopac_example_n1'])
    ens, poss = read_report('all_equations_mnw_3a4_rr', xyz=True)
    # ens, poss = read_report('equations_mopac_example_n0', xyz=True)
    to_names = poss.pop(0)
    names = [i[0] for i in to_names]
    atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
    interaction = length_xyz_interaction((names, atoms))
    ss = system.solve()
    mops = [ens[0]]
    appr = [apply_solution(ss, interaction)]
    for en, item in zip(ens, poss):
        atoms = [np.array(list(map(float, i[1::]))) for i in item]
        interaction = length_xyz_interaction((names, atoms))
        if isinstance(interaction, dict) and system.can_i_solve_it(interaction):
            appr.append(apply_solution(ss, interaction))
            mops.append(en)

    appr = np.array(appr)
    mops = np.array(mops)
    # divider = Spherical_divider(n)
    # genetic_to_the_aim(Molecule(file_name, n, divider), Molecule(to_file, n, divider), write=True,
    #  file_log=saver, trash_keeper=trasher)

    # print(np.mean([abs(i-j) for i, j in zip(mops, appr)]))
    print(np.corrcoef(np.array(mops), np.array(appr)))
    plt.xlabel('Mopac energy, eV')
    plt.ylabel('Approximate energy, eV')
    plt.scatter(mops, appr)
    plt.show()

    # c = 1 if system.find_equation_in_system(interaction) != -1 else 0
    # for en, item in zip(ens, poss):
    #     atoms = [np.array(list(map(float, i[1::]))) for i in item]
    #     interaction = length_xyz_interaction((names, atoms))
    #     c += 1 if system.find_equation_in_system(interaction) != -1 else 0
    # print(c)


    # # file_name = './prepared_mols2/3a_opted.mol2'
    # # to_file = './prepared_mols2/4_opted.mol2'
    # file_name = './ordered_mol2/js_exapmle_init.mol2'
    # to_file = './ordered_mol2/js_exapmle_finish.mol2'
    # saver = 'equations_mopac_example2'
    # # keep_equations_path(file_name, to_file, saver)
    # read_eqs_path(file_name, to_file, saver, n=20, solve_rel_mean=False)

    def procedure_for_correlate_energies():
        # file_name = './ordered_mol2/js_exapmle_init.mol2'
        file_name = './prepared_mols2/3a_opted.mol2'
        # to_file = './ordered_mol2/js_exapmle_finish.mol2'
        to_file = './prepared_mols2/4_opted.mol2'
        saver = 'equations_mnw_3a4_4'
        trasher = 'all_eqs_0'
        n = 20
        system = Equation_system()
        # system.from_reports(['equations_mn_3a4_0', 'equations_mn_3a4_1'])
        system.from_reports(['equations_mopac_example_n0', 'equations_mopac_example_n1'])
        # ens, poss = read_report('equations_mn_3a4_2', xyz=True)
        ens, poss = read_report('equations_mopac_example_n0', xyz=True)
        to_names = poss.pop(0)
        names = [i[0] for i in to_names]
        atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
        interaction = length_xyz_interaction((names, atoms))
        ss = system.solve()
        mops = [ens[0]]
        appr = [apply_solution(ss, interaction)]
        for en, item in zip(ens, poss):
            atoms = [np.array(list(map(float, i[1::]))) for i in item]
            interaction = length_xyz_interaction((names, atoms))
            if isinstance(interaction, dict) and system.can_i_solve_it(interaction):
                appr.append(apply_solution(ss, interaction))
                mops.append(en)

        appr = np.array(appr)
        mops = np.array(mops)
        # it's for new pathes
        # divider = Spherical_divider(n)
        # genetic_to_the_aim(Molecule(file_name, n, divider), Molecule(to_file, n, divider), write=True,
        #  file_log=saver, trash_keeper=trasher)

        print(np.mean([abs(i - j) for i, j in zip(mops, appr)]))
        print(np.corrcoef(np.array(mops), np.array(appr)))
        plt.xlabel('Mopac energy, eV')
        plt.ylabel('Approximate energy, eV')
        plt.scatter(mops, appr)
        plt.show()

