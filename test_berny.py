from pytest import approx
import numpy as np
from berny import optimize, geomlib
from berny.solvers import MopacSolver

ethanol = geomlib.loads("""\
9
Ethanol
  H      1.8853     -0.0401      1.0854
  C      1.2699     -0.0477      0.1772
  H      1.5840      0.8007     -0.4449
  H      1.5089     -0.9636     -0.3791
  C     -0.2033      0.0282      0.5345
  H     -0.4993     -0.8287      1.1714
  H     -0.4235      0.9513      1.1064
  O     -0.9394      0.0157     -0.6674
  H     -1.8540      0.0626     -0.4252
""", 'xyz')


def test_basic():
    solver = MopacSolver()
    final = optimize(solver, ethanol, steprms=0.01, stepmax=0.05)
    inertia_princpl = np.linalg.eigvalsh(final.inertia)
    print('eigvalsh', inertia_princpl)

# test_basic()