import numpy as np
import matplotlib.pyplot as plt

def mopac_gradien_plotter(file_name): #'5_opt.out'
    grads = []
    with open(file_name, 'r') as f:
        next(l for l in f if 'Geometry optimization using L-BFGS' in l)
        line = f.readline()
        while len(line.split())>0:
            # print(line)
            grads.append(float(line.split()[-3]))
            line = f.readline()

    plt.plot(list(range(len(grads))), grads)
    plt.xlabel('Num of step')
    plt.ylabel('Step gradient')
    plt.show()

def log_berny_plotter(file_name):
    log = open(file_name, 'r').read()
    log = log.split('\n')
    loge = list(filter(lambda x: 'Energy:' in x, log))
    logens = [float(i.split()[-1]) for i in loge]
    plt.plot(list(range(len(logens))), logens)
    plt.xlabel('Num of step')
    plt.ylabel('Step energy')
    plt.show()

if __name__ == '__main__':
    # log_berny_plotter('/home/anastasiia/PycharmProjects/chem_structure/5_steprms_0.01_stepmax_0.05.log')
    log_berny_plotter('/home/anastasiia/PycharmProjects/chem_structure/5_steprms_0.05_stepmax_0.1.log')