import numpy as np
import random
import matplotlib.pyplot as plt
from mopac_worker import get_energy_of_xyz, mopacOut_to_xyz_with_energy,\
    writeInputFileFromXYZ
import shutil, tempfile
from subprocess import call
from os import path
# from collections import Counter
alias = '/opt/mopac/run_script.sh'


def spoil(points, d=0.05, diffs=False):
    spoileds = []
    ds = []
    for i in points:
        adv = np.array([random.normalvariate(0, d)
                                  for _ in range(3)])
        spoileds.append(np.array(i) + adv)
        if diffs:
            ds.append(np.linalg.norm(adv))
    if diffs: return spoileds, sum(ds)
    return spoileds


def spoil_de_dg(file_name, d=0.5, n=250):
    tmpdir = tempfile.mkdtemp()
    lines = open(file_name, 'r').read().split('\n')
    points = [np.array([float(i) for i in line.split('\t')[1::]]) for line in lines[2::]]

    original_energy = round(get_energy_of_xyz(file_name), 2)
    for _ in range(n):
        spoileds, diff = spoil(points, d=d, diffs=True)
        cur_f = path.join(tmpdir, file_name[:-4:]+'_spoiled.xyz')
        with open(cur_f, 'w') as f_spoiled:
            f_spoiled.write(str(len(points)) + '\n')
            for i in range(len(points)):
                f_spoiled.write('\n{}\t{}\t{}\t{}'.format(lines[i + 2].split('\t')[0], *spoileds[i]))
        options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Single Point',
                   'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}

        ens = round(get_energy_of_xyz(cur_f), 2)
        plt.scatter(diff, ens - original_energy, color='b')
    plt.title('With opt ' + file_name[:-4:]+'_state std=' + str(d))
    plt.xlabel('Geoemtrical diffrence, A')
    plt.ylabel('Energy diffence')
    plt.show()
    shutil.rmtree(tmpdir)


def spoil_and_opt(file_name, d=0.5, n=250):
    tmpdir = tempfile.mkdtemp()
    lines = open(file_name, 'r').read().split('\n')
    points = [np.array([float(i) for i in line.split('\t')[1::]]) for line in lines[2::]]

    original_energy = round(get_energy_of_xyz(file_name), 2)
    for _ in range(n):
        spoileds, diff = spoil(points, d=d, diffs=True)
        cur_f = path.join(tmpdir, file_name[:-4:]+'_spoiled.xyz')
        with open(cur_f, 'w') as f_spoiled:
            f_spoiled.write(str(len(points)) + '\n')
            for i in range(len(points)):
                f_spoiled.write('\n{}\t{}\t{}\t{}'.format(lines[i + 2].split('\t')[0], *spoileds[i]))
        options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
                   'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}
        opted = cur_f[:-4:]+'_opt'
        writeInputFileFromXYZ(options, cur_f,  opted+'.mop')

        call([alias, opted+'.mop'])
        mopacOut_to_xyz_with_energy(opted, opted+'_opted.xyz')
        ens = round(get_energy_of_xyz(opted+'_opted.xyz'), 2)
        # print(get_energy_of_xyz(cur_f))
        # print(ens)
        plt.scatter(diff, ens - original_energy, color='b')
    plt.title(file_name[:-4:]+'_state std=' + str(d))
    plt.xlabel('Geoemtrical diffrence, A')
    plt.ylabel('Energy diffence')
    plt.show()
    shutil.rmtree(tmpdir)


if __name__ == '__main__':
    # spoil_and_opt('4a_opt.xyz', d=0.05, n=100)
    spoil_and_opt('5_opt.xyz', d=0.05, n=100)
    # spoil_de_dg('5_opt.xyz', d=0.05, n=100)
    # print(get_energy_of_xyz('4a_opt.xyz'))
    # print(get_energy_of_xyz('5_opt.xyz'))
