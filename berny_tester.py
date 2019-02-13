from berny import Berny, optimize, geomlib
from berny.solvers import MopacSolver
import numpy as np
# import sys

# orig_stdout = sys.stdout
srm, smax = 0.05, 0.1
fxyz, logf, fopt = 'TS-5-6-19.xyz', 'TS-5-6_steprms_'+str(srm)+'_stepmax_'+str(smax)+'.log', '5_steprms_'+str(srm)+'_stepmax_'+str(smax)+'.xyz'
# f = open(logf, 'w')
# sys.stdout = f

optimizer = Berny(geomlib.readfile(fxyz), steprms=srm, stepmax=smax, maxsteps=150)
final = optimize(optimizer, MopacSolver(cmd='/opt/mopac/run_script.sh'))
inertia_princpl = np.linalg.eigvalsh(final.inertia)
final.dump(open('TS-5-6_bopt.xyz', 'w'), 'xyz')

#f.close()

# for geom in optimizer:
    # get energy and gradients for geom
    # optimizer.send((energy, gradients))
