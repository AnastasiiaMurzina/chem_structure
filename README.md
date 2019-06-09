This one allows show many possible compounds.

Compounds are different if at least one bonds has other energy. Molecule can be kept with Molecule('*.mol2')
 and also it probably used with divider:
###############################
```from mol_api import Molecule
from quadro_with_rotate_class import Spherical_divider

n = 15
divider = Spherical_divider(n=n) #choose accuracy: recommended more than 12

reactant = Molecule('my_reactant.mol2', n=n, divider=divider)
prod = Molecule('my_product.mol2', n=n, divider=divider) #choose the same n
```
####################### 

Now all functions of getting energy are connected with MOPAC[1], but it'll be good if it'll possible to
use another force field calculation programs.


The reaction path could be probed with two different released methods:
-- if the product is known:

```
from mixed_searcher import approx_genetic_to_the_aim, genetic_to_the_aim

genetic_to_the_aim(reactant, prod, file_log=file_log) #energies'll be logged in file
```
-- and not known (series of mutation and than it could be optimized):
```
reactant.mutation() # 3 types of mutation: add/del bond, change length of bond, change section of bond
```
### Many molecules -- SOON
Now it've released like '0'-bonds in mol2-file

### Geometric align

Using rmsd (pip install rmsd) 
```
import rmsd

coord -= rmsd.centroid(coord)
coord3 -= rmsd.centroid(coord3)
rotate = rmsd.kabsch(coord3, coord)
coord3 = np.dot(coord3, rotate)
rmsd.rmsd(coord, coord3)
```

1. Stewart, J. J. (1990). MOPAC: a semiempirical molecular orbital program. Journal of computer-aided molecular design, 4(1), 1-103.