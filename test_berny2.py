from berny.solvers import MopacSolver
from berny2.berny import Berny, geomlib
from berny2.berny.solvers import MopacSolver

optimizer = Berny(geomlib.readfile('mols_dir/Ethene.xyz'), debug=True)
solver = MopacSolver()
for geom in optimizer:
    energy, gradients = solver(geom)
    optimizer.send((energy, gradients))
