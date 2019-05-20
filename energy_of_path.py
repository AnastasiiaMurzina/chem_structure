from os import path
from tempfile import mkdtemp
from shutil import rmtree

import matplotlib.pyplot as plt
from mopac_worker import get_heat_of_xyz, get_energy_of_xyz

def read_report(filename):
    with open(filename, 'r') as f:
        n_steps = int(f.readline())
        ll = []
        for _ in range(n_steps):
            line = f.readline()
            if line == '':
                break
            n_mols = int(line)
            f.readline()
            lines = [f.readline() for _ in range(n_mols)]
            f.readline()
            ll.append(lines)
    return ll

def get_func_of_report(function, filename):
    ll = read_report(filename)
    tmp = mkdtemp()
    xyzf = path.join(tmp, 'xyzf.xyz')
    ens = []
    for positions in ll:
        with open(xyzf, 'w') as f:
            f.write(str(len(positions))+'\n\n')
            for position in positions:
                f.write(position)
        ens.append(function(xyzf, tmpdir=tmp))
    rmtree(tmp)
    return ens


if __name__ == '__main__':
    from mol_api import Molecule, random_to_the_aim_search
    n=15
    for i in range(16, 30):
        ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n)
        pr = Molecule('./prepared_mols2/4_opted.mol2', n=n)
        pr.refresh_dimensional()
        report_name = '3a->4_report_'+str(i)
        ms = random_to_the_aim_search(ln, pr, write=True, file_log= report_name)
        h = get_func_of_report(get_energy_of_xyz, report_name)
        plt.title('3a->4 mopac energies, '+report_name)
        print(h)
        print(max(h))
        plt.plot(list(range(len(h))), h)
        plt.xlabel('number of step')
        plt.ylabel('kcal/mol')
        plt.show()
