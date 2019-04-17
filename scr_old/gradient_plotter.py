import numpy as np
import matplotlib.pyplot as plt
import tempfile, shutil, subprocess

def mopac_gradien_plotter(file_name): #'5_opt.out'
    grads = []
    with open(file_name, 'r') as f:
        next(l for l in f if 'Geometry optimization using L-BFGS' in l)
        line = f.readline()
        while len(line.split()) > 0:
            # print(line)
            grads.append(float(line.split()[-3]))
            line = f.readline()

    plt.plot(list(range(len(grads))), grads)
    plt.xlabel('Num of step')
    plt.ylabel('Step gradient')
    plt.show()

def mopac_aux_energy_getter(aux):
    energies = []
    with open(aux.split('.')[-2]+'.mop', 'r') as fnames:
        text = fnames.read().split('\n')
        text = list(filter(lambda x: x!='', text))
        names = [i.split()[0] for i in text[2::]]
    td = tempfile.mkdtemp()
    with open(aux, 'r') as f:
        i = 1
        line = ''
        while not('#' in line):
            next(l for l in f if '_UPDATED:ANGSTROMS' in l)
            xyz = []
            for _ in range(len(names)):
                line = f.readline().split()
                xyz.append(line)
            to_calc_en = td+'/'+str(i)+'.mop'
            with open(to_calc_en, 'w') as calc:
                calc.write(' AUX LARGE CHARGE=0 SINGLET NOOPT PM6\nTitle\n\n')
                for ix in zip(names, xyz):
                    calc.write('{}\t{}\t0\t{}\t0\t{}\t0\n'.format(ix[0], *ix[1]))
            subprocess.call(['/opt/mopac/run_script.sh', to_calc_en])
            with open(td+'/'+str(i)+'.out', 'r') as fen:
                line = next(l for l in fen if 'TOTAL ENERGY' in l)
                energies.append(float(line.split()[-2]))
            i += 1
            line = next(f)
    shutil.rmtree(td)
    return energies

def mopac_aux_energy_plotter(aux):
    ens = mopac_aux_energy_getter(aux)
    plt.plot(list(range(len(ens))), ens)
    plt.xlabel('Num of step')
    plt.ylabel('Step energy')
    plt.title('Mopac energy optimizer')
    plt.show()


def log_berny_plotter(file_name, mopac_energy=True, title='Berny energy optimizer'):
    log = open(file_name, 'r').read()
    log = log.split('\n')
    loge = list(filter(lambda x: 'Energy:' in x, log))
    if mopac_energy:
        logens = [float(i.split()[-1])*27.21139 for i in loge] #hartree to Ev
    else:
        logens = [float(i.split()[-1]) for i in loge]
    plt.plot(list(range(len(logens))), logens)
    plt.xlabel('Num of step')
    plt.ylabel('Step energy')
    plt.title(title)
    plt.ticklabel_format(useOffset=False)
    plt.show()

if __name__ == '__main__':
    # mopac_aux_energy_plotter('TS-5-6_ts.aux')
    # log_berny_plotter('/home/anastasiia/PycharmProjects/chem_structure/5_steprms_0.01_stepmax_0.05.log')
    # log_berny_plotter('6_bopt.log')
    pass