import matplotlib.pyplot as plt
import copy
import numpy as np
from tempfile import mkdtemp
from shutil import rmtree

from mol_api import Molecule, compare_structers
from layouter import dimensional_structure


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
        compared = pr.compare_with(np.array([i for _, i in ch_dim.items()]))
        _, s_comp_c, _ = st.notation.diff(pr.notation)
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

def genetic_to_the_aim(reactant, product, write=False,
                             file_log='reaction_report'):
    """
    :return: energies_path and length of sections changes
    """
    def apply_change():
        mutant.refresh_dimensional()
        en = mutant.get_energy()
        if en is None:
            return False
        delta = en - path[-1]
        if np.random.random() < np.exp(-delta/10.):
            path.append(path[-1]+delta)
            if write:
                mutant.to_xyz(file_log, title=str(path[-1]), mode='a')
                with open(file_log, 'a') as f_w:
                    f_w.write('################\n')
            return True
        return False

    with open(file_log, 'w') as f_w:
        f_w.write(str(len(reactant.notation.s_diff(product.notation))
                      +len(reactant.notation.l_diff(product.notation)))+'\n')
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
    while reactant.notation.l_change_step(product.notation) != -1:
        mutant = copy.deepcopy(reactant)
        if apply_change(): reactant = copy.deepcopy(mutant)
    while reactant.notation.s_change_step(product.notation) != -1:
        mutant = copy.deepcopy(reactant)
        if apply_change(): reactant = copy.deepcopy(mutant)
    return path #, msds

    #########two halves of random path#########
    # s = reactant.notation.s_change(product.notation, follow_energy=True, sh=True)
    # l = reactant.notation.l_change(product.notation, follow_energy=True, sh=True)
    # return s+l, len(s)
    ###########################################


class Equation_system:
    def __init__(self, first_linear, energy):
        self.equations = [first_linear]
        self.variables = set(first_linear.keys())
        self.energy = [energy]

    def push(self, linear, energy):
        self.equations.append(linear)
        self.variables.update(linear.keys())
        self.energy.append(energy)

    def check_system(self):
        return len(self.energy) > len(self.variables)

    def can_i_solve_it(self, variables):
        return len(variables.keys() - self.variables) == 0

    def to_matrix(self):
        x = np.zeros((len(self.energy), len(self.variables)))
        order_of_vars = []
        for ix_e, var in enumerate(self.variables):
            order_of_vars.append(var)
            for ix, i in enumerate(self.equations):
                x[ix][ix_e] = i.get(var, 0)
        return x, np.array(self.energy), order_of_vars

    def solve(self):
        system, ens, order = self.to_matrix()
        print(system, ens)
        try:
            print('solution', np.linalg.lstsq(system, ens))
            return {i: j for i, j in zip(order, np.linalg.lstsq(system, ens))}
        except np.linalg.LinAlgError:
            return -1

def apply_solution(solution: dict, interacted: dict):
    print(solution)
    print(len(solution))
    en = 0
    for key, item in interacted.items():
        en += solution[1][key]*item
    return en



def multi_criterical(mol1: Molecule, to_mol: Molecule, n=1000, write=False, file_log='multi_report', limits=2):
    temp_dir = mkdtemp()

    def apply_change():
        mut = copy.deepcopy(mut_pr)
        mut.refresh_dimensional()

        distance.append(msd_pr)
        difference.append(difference_pr)
        energies.append(energy_pr)
        if write:
            mut.to_xyz(file_log, title=str(energy_pr), mode='a')
            with open(file_log, 'a') as f_w:
                f_w.write('################\n')

    with open(file_log, 'w') as f_w:
        f_w.write(str(n)+'\n')

    mut = copy.deepcopy(mol1)
    distance = [compare_structers(mol1.to_positions_array(), to_mol.to_positions_array())]
    difference = [mol1.notation.diff(to_mol.notation)]
    energies = [mol1.get_energy(tmp=temp_dir)]
    solver = Equation_system(mut.interact_pair(), energies[0])
    print('initial parameters')
    print('initial energy', energies[0])
    print('initial rmsd', distance[0])
    print('initial differnce', difference[0])

    for _ in range(n):
        mut_pr = mut.mutation(probailities=[0.5, 1, 1])
        mut_pr.refresh_dimensional()

        difference_pr = mut_pr.notation.diff(to_mol.notation) # d[1] - diff_of_sections, d[2] - summ diff of lengths
        msd_pr = compare_structers(mut_pr.to_positions_array(), to_mol.to_positions_array())
        vars_linear_approx = mut_pr.interact_pair()
        # if np.random.random() < np.exp(-len(solver.variables)/len(solver.equations))\
        #         and solver.check_system() and solver.can_i_solve_it(vars_linear_approx):
        #     print('could')
        #     ss = solver.solve()
        #     if ss == -1:
        #         energy_pr = None
        #     else:
        #         energy_pr = apply_solution(ss, vars_linear_approx)
        # else:
        energy_pr = mut_pr.get_energy(tmp=temp_dir)
        if energy_pr is None:
            continue
        # solver.push(mut_pr.interact_pair(), energy_pr)
        # if difference[-1] < difference_pr or distance[-1] < msd_pr\
        #         or np.random.random() < np.exp(-(energy_pr/energies[-1])**2):
        if difference[-1] < difference_pr or distance[-1] < msd_pr \
                and energy_pr < -4680.:
            apply_change()
            print('rmsd accept', msd_pr)
            print('difference accept', difference_pr)
            print('energy accept', energy_pr)

    print(difference)
    print(distance)
    print(energies)


if __name__ == '__main__':
    n = 2
    # reaction = 'mopac_example' # '3a->4' #
    reaction = '3a->4' #
    # reaction = 'vanadii'
    if reaction == '3a->4':
        ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n)
        pr = Molecule('./prepared_mols2/4_opted.mol2', n=n)
        pr.refresh_dimensional()
    # elif reaction == 'vanadii':
    #     ln = Molecule('./vanadii/3a_singlet_opted.mol2', n=n)
    #     pr = Molecule('./vanadii/ts_3a_4a_opted.mol2', n=n)
    #     pr.refresh_dimensional()
    # else:
    #     ln = Molecule('./ordered_mol2/js_exapmle_init.mol2', n=n)
    #     pr = Molecule('./ordered_mol2/js_exapmle_finish.mol2', n=n)
    #     pr.refresh_dimensional()
    lstl_solver = Equation_system(ln.interact_pair(), ln.get_energy())
    multi_criterical(ln, pr, n=1000)

    # ln.notation.s_change_step()
    # maxs = []
    # for i in range(5):
    #     ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n)
    #     pr = Molecule('./prepared_mols2/4_opted.mol2', n=n)
    #     pr.refresh_dimensional()
    #     ms = genetic_to_the_aim(ln, pr, write=True, file_log=reaction+'_gen_report_s20_'+str(i))
    #     # print(ms, )
    #     maxs.append(max(ms))
    # print(maxs)
    # print(min(maxs))

    # print(ln.notation.get_energy())
    # print(pr.notation.get_energy())
    # real_random_path(ln, pr, n=10000, write=True, file_log=reaction+'_random_report')
    # print(max(p))
    # print(p)
    # print(max(ms))
    # print(ms)
    # print(ln.notation.l_change(pr.notation, follow_energy=True))
    # print(ln.notation.s_change(pr.notation, follow_energy=True))


    # write_mol2_file('My_' + name + '_' + 'q0.mol2', atoms_info, dim_structure, bonds=paired)