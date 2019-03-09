from berny import Berny, optimize, geomlib
from berny.solvers import MopacSolver
from gradient_plotter import log_berny_plotter
import numpy as np
import sys

# orig_stdout = sys.stdout
srm, smax = 0.05, 0.1
fxyz, logf, fopt = '/home/anastasiia/PycharmProjects/chem_structure/opted_mols/3-MnH2_opted.xyz',\
                   '3_opted'+str(srm)+'_stepmax_'+str(smax)+'.log', '3_opted_'+str(srm)+'_stepmax_'+str(smax)+'.xyz'
f = open(logf, 'w')
sys.stdout = f

optimizer = Berny(geomlib.readfile(fxyz), steprms=srm, stepmax=smax, maxsteps=150)
final = optimize(optimizer, MopacSolver(cmd='/opt/mopac/run_script.sh'))
inertia_princpl = np.linalg.eigvalsh(final.inertia)
final.dump(open('3_berny_opt_mopac.xyz', 'w'), 'xyz')
# final.dump(open('TS-5-6_bopt.xyz', 'w'), 'xyz')
f.close()
log_berny_plotter(logf, title='3_berny_opt_mopac')

# for geom in optimizer:
    # get energy and gradients for geom
    # optimizer.send((energy, gradients))
