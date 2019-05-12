It's the structure for solving the stereochemistry combinatoric problem. This one allows show all possible stereo isomers.

Compouds are different if at least one bonds has other energy.

The reason may be another order of elements (if there'are 4 bonds: append clockwise and counter-clockwise), may be angles, etc.

Isomer must be a stable compoud.
###############################

from berny import Berny, geomlib

optimizer = Berny(geomlib.readfile('mols_dir/Aniline.xyz'))
for geom in optimizer:
    # get energy and gradients for geom
    optimizer.send((energy, gradients))


####################### 

python -mSimpleHTTPServer 8123



### One molecule from original mol2 to restored from notation

```
from mol_api import Molecule
mol = Molecule("example.mol2")
```

### Many molecules -- SOON
Now it've released like '0'-bonds in mol2-file
### Align

Using rmsd (pip install rmsd) 
```
import rmsd

coord -= rmsd.centroid(coord)
coord3 -= rmsd.centroid(coord3)
rotate = rmsd.kabsch(coord3, coord)
coord3 = np.dot(coord3, rotate)
rmsd.rmsd(coord, coord3)
```