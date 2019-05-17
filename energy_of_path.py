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
    # heats = get_heat_of_report('reaction_report')
    # # heats = [465.80044, 465.80044, 465.80044, 463.09968, 463.09968, 463.09968, 463.09968, 442.92165, 458.45804, 458.45804, 458.45804, 458.45804, 458.45804, 458.45804, 458.45804, 478.07156, 521.10831, 521.10831, 521.10831, 521.10831, 521.10831, 529.20735, 527.3009, 527.3009, 527.3009, 527.3009, 546.26281, 546.26281, 609.47926, 609.47926, 610.5946, 615.73759, 573.78185, 480.75425, 442.52887, 455.42596]
    # plt.plot(list(range(len(heats))), heats)
    # plt.xlabel('number of step')
    # plt.ylabel('kcal/mol')
    # plt.show()
    h = get_func_of_report(get_energy_of_xyz, '3a->4_report_opted')
    plt.title('3a->4 mopac heats')
    # h = get_heat_of_report('mopac_example_report')
    print(h)
    print(max(h))
    plt.plot(list(range(len(h))), h)
    plt.xlabel('number of step')
    plt.ylabel('kcal/mol')
    plt.show()