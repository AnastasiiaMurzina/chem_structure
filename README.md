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
