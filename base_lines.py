from os import path
import mol2_worker
from mopac_worker import get_energy_of_xyz, get_heat_of_xyz
import tempfile, shutil
from thermo_characteristics import get_dG
import numpy as np
import matplotlib.pyplot as plt
import rmsd

def get_linear_path_energies(reactant, names, product, n=10):
    vbegin = [np.array(f) for f in reactant]
    vend = [np.array(f) for f in product]
    vector_to_product = [f - t for f, t in zip(vend, vbegin)]
    mols = len(vbegin)
    x, y = list(range(n+1)), []
    tmp = tempfile.mkdtemp()
    tmp_file = path.join(tmp, 'current.xyz')
    for ind in x:
        current_state = [f + vp * ind / n for f, vp in zip(vbegin, vector_to_product)]
        with open(tmp_file, 'w') as f:

            f.write(str(mols) + '\n\n')
            for inx in range(mols):
                f.write('{}\t{}\t{}\t{}\n'.format(names[inx], *current_state[inx]))
        y.append(get_energy_of_xyz(tmp_file, tmpdir=tmp))
    shutil.rmtree(tmp)
    return x, y

def get_linear_path_heats(reactant, product, n=10):
    vbegin = [np.array(f[1::]) for f in reactant]
    vend = [np.array(f[1::]) for f in product]
    vector_to_product = [f - t for f, t in zip(vend, vbegin)]
    mols = len(vbegin)
    x, y = list(range(n+1)), []
    tmp = tempfile.mkdtemp()
    tmp_file = path.join(tmp, 'current.xyz')
    for ind in x:
        current_state = [f + vp * ind / n for f, vp in zip(vbegin, vector_to_product)]
        with open(tmp_file, 'w') as f:
            f.write(str(mols) + '\n\n')
            for inx in range(mols):
                f.write('{}\t{}\t{}\t{}\n'.format(reactant[inx][0], *current_state[inx]))
        y.append(get_heat_of_xyz(tmp_file, tmpdir=tmp))
    shutil.rmtree(tmp)
    return x, y

def get_linear_path_gibbs_ezpe(reactant, product, n=10):
    vbegin = [np.array(f[1::]) for f in reactant]
    vend = [np.array(f[1::]) for f in product]
    vector_to_product = [f - t for f, t in zip(vend, vbegin)]
    mols = len(vbegin)
    x, g, zpe = list(range(n+1)), [], []
    tmp = tempfile.mkdtemp()
    tmp_file = path.join(tmp, 'current.xyz')
    for ind in x:
        current_state = [f + vp * ind / n for f, vp in zip(vbegin, vector_to_product)]
        with open(tmp_file, 'w') as f:
            f.write(str(mols) + '\n\n')
            for inx in range(mols):
                f.write('{}\t{}\t{}\t{}\n'.format(reactant[inx][0], *current_state[inx]))
        _, _, gibbs, ez = get_dG(tmp_file, tmp=tmp)
        g.append(gibbs/1000.)
        zpe.append(ez)
    shutil.rmtree(tmp)
    return x, g, zpe

def show_path_energies(steps, energies, product_name='', reactant_name='', article_TS_energy=[]):
    plt.plot(steps, energies)
    if article_TS_energy != []:
        plt.plot(steps, article_TS_energy)
    plt.title(reactant_name+'->'+product_name+' linear path in energies')
    plt.ylabel('EV')
    plt.show()


def linear_line(reactant, product, product_name='', reactant_name=''):
    vbegin = [np.array(f[1::]) for f in reactant]
    vend = [np.array(f[1::]) for f in product]
    vector_to_product = [f - t for f, t in zip(vend, vbegin)]
    n = 10
    mols = len(vbegin)
    x, y = [], []
    tmp = tempfile.mkdtemp()
    tmp_file = path.join(tmp, 'current.xyz')
    for ind in range(n):
        x.append(ind)
        current_state = [f + vp * ind / n for f, vp in zip(vbegin, vector_to_product)]
        with open(tmp_file, 'w') as f:
            f.write(str(mols) + '\n\n')
            for inx in range(mols):
                f.write('{}\t{}\t{}\t{}\n'.format(reactant[inx][0], *current_state[inx]))
        _, _, gibbs, zpe = get_dG(tmp_file, tmp=tmp)
        y.append(gibbs/1000.-zpe)
    plt.plot(x, y)
    # energy = get_energy_of_xyz('/home/anastasiia/PycharmProjects/chem_structure/article_xyz/TS-3a-4.xyz')
    # plt.plot(x, [energy]*len(x))
    plt.xlabel('steps')
    plt.ylabel('kcal/mol')
    plt.title(reactant_name+'->'+product_name+' linear path in free Gibbs energy')
    plt.show()
    print('free Gibbs energy - E_zpe threshold of reaction ', max(y) - y[0])
    print('free Gibbs energy - E_zpe delta of reaction', y[-1] - y[0])
    shutil.rmtree(tmp)

    def mc_rmsd_line(reactant, product):
        vthree = [np.array(f[1::]) for f in reactant]
        vfour = [np.array(f[1::]) for f in product]

        n = 10
        mols = 51
        x, y = [], []
        tmp = tempfile.mkdtemp()
        tmp_file = path.join(tmp, 'current.xyz')

        for ind in range(n):
            x.append(ind)
            # current_state = [f + vp * ind / n for f, vp in zip(vthree, vector_to_product)]
            current_state = [] #TODO get notation, change notation, accept or reject changes
            with open(tmp_file, 'w') as f:
                f.write(str(mols) + '\n\n')
                for inx in range(mols):
                    f.write('{}\t{}\t{}\t{}\n'.format(reactant[inx][0], *current_state[inx]))
            y.append(get_energy_of_xyz(tmp_file))
        # scaled = 103*(y[-1] - y[0])
        # y =[yy*scaled for yy in y]
        plt.plot(x, y)
        plt.title('3a->4 MC_RMSD path in Jol')
        plt.show()

        shutil.rmtree(tmp)



if __name__ == "__main__":

    # a = mol2_worker.Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/3a_opted.mol2')
    # three_a = a.to_positions_array()
    # four_mol = mol2_worker.Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/4_opted.mol2')
    # four = four_mol.to_positions_array()




    three_a, n = mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/3a_opted.xyz', names=True)
    four = mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/4_opted.xyz')
    ts_3a_4 =  mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/article_xyz/TS-3a-4.xyz')

    # three_a, n = mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/js_exapmle_init.xyz', names=True)
    # four = mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/js_exapmle_finish.xyz')
    # ts_3a_4 =  mol2_worker.xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/article_xyz/TS-3a-4.xyz')

    # seven = mol2_worker.xyz_to_array('./ordered_mol2/7_opted.xyz')
    # eight = mol2_worker.xyz_to_array('./ordered_mol2/8-Mn-ads-EtOH_opted.xyz')
    # linear_line(three_a, four, reactant_name='3a', product_name='4')
    # linear_line(three_a, four, reactant_name='3a', product_name='4')
    x, y = get_linear_path_energies(three_a, n, four, n=515)
    print(x)
    print(y)
    # show_path_energies(list(range(y)), y, reactant_name='3a', product_name='4', article_TS_energy=[-4704.37]*30)

    # linear_line(seven, eight, reactant_name='7', product_name='8')


