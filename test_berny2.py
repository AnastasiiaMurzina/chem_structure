from berny.solvers import MopacSolver
import numpy as np
from berny2.berny import Berny, geomlib, optimize
from berny2.berny.solvers import MopacSolver

mol = geomlib.readfile('mols_dir/Ethanol.xyz')
solver = MopacSolver()
final = optimize(solver, mol, steprms=0.01, stepmax=(1+np.cos(0.4*np.pi))**0.5)
inertia_princpl = np.linalg.eigvalsh(final.inertia)
print(final.coords)
# print(solver)
# for i in solver:
    # print(i)
# for geom in optimizer:
#     print(geom)
    # energy, gradients = solver(geom)
    # optimizer.send((energy, gradients))
