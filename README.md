It's the structure for solving the stereochemistry combinatoric problem. This one allows show all possible stereo isomers.

Compouds are different if at least one bonds has other energy.

The reason may be another order of elements (if there'are 4 bonds: append clockwise and counter-clockwise), may be angles, etc.

Isomer must be a stable compoud.


### One molecule from original mol2 to restored from notation

```
from mol2_worker import atoms_and_bonds, xyz_names_bonds  
from mol2_chain import mol2_to_notation, dimensional_structure, write_mol2_file

atoms_info = atoms_and_bonds(name + '.mol2') #get information about atoms from mol2 - it's needed for write new mol2 later
notation_keeper = mol2_to_notation(xyz_names_bonds(name + '.mol2'), method='ten') #get notation
paired = bonds_of_paired(notation_keeper[1])
dim_structure = dimensional_structure([notation_keeper[0], paired], method='ten') #return to dimentional xyz
write_mol2_file('new_name+'.mol2', atoms_info, dim_structure, bonds=paired) #write restored dimentional structure
```

### Many molecules
```
from mol2_worker import atoms_and_bonds, xyz_names_bonds 
from many_mols import get_notation_many_mols, places_for_zero_bonds, zero_connecter, unpack_with_zero_bonds, write_mol2_file

bs, ass = xyz_names_bonds(name + '.mol2')
atoms_notation, bonds_notation = get_notation_many_mols(ass, bs)

needs_to_zero_discribe = places_for_zero_bonds(ass, bs)

zero_bonds = zero_connecter(ass, needs_to_zero_discribe)
coords = unpack_with_zero_bonds(atoms_notation, bonds_notation, zero_bonds)

write_mol2_file('new_file_name.mol2', ass, coords, to_two_ways_bond(bs, with_attr=True))
```
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



