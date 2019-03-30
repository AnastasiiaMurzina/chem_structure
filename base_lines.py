from os import path
import mol2_worker
from mopac_worker import get_energy_of_xyz
import tempfile, shutil
from thermo_characteristics import get_dG
import numpy as np
import matplotlib.pyplot as plt
import rmsd

if __name__ == "__main__":

    a = mol2_worker.Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/3a_opted.mol2')
    # three_a = a.to_positions_array()
    four_mol = mol2_worker.Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/4_opted.mol2')
    # four = four_mol.to_positions_array()


    three_a = xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/3a_scaled_opted.xyz')
    four = xyz_to_array('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/4_opted.xyz')

    def linear_line(reactant, product):
        vthree = [np.array(f[1::]) for f in reactant]
        vfour = [np.array(f[1::]) for f in product]
        vector_to_product = [f-t for f, t in zip(vfour, vthree)]
        n = 10
        mols = 51
        x, y = [], []
        tmp = tempfile.mkdtemp()
        tmp_file = path.join(tmp, 'current.xyz')
        for ind in range(n):
            x.append(ind)
            current_state = [f + vp * ind / n for f, vp in zip(vthree, vector_to_product)]
            with open(tmp_file, 'w') as f:
                f.write(str(mols)+'\n\n')
                for inx in range(mols):
                    f.write('{}\t{}\t{}\t{}\n'.format(reactant[inx][0], *current_state[inx]))
            y.append(get_energy_of_xyz(tmp_file)/27.2114*2625.5)
        plt.plot(x, y)
        plt.title('3a->4 linear path in Jol')
        plt.show()
        print(y[-1]-y[0])
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
            y.append(get_energy_of_xyz(tmp_file) / 27.2114 * 2625.5)
        plt.plot(x, y)
        plt.title('3a->4 MC_RMSD path in Jol')
        plt.show()
        print(y[-1] - y[0])
        shutil.rmtree(tmp)



