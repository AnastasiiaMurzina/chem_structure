It's the structure for solving the stereochemistry combinatoric problem. This one allows show all possible stereo isomers.

Compouds are different if at least one bonds has other energy.

The reason may be another order of elements (if there'are 4 bonds: append clockwise and counter-clockwise), may be angles, etc.

Isomer must be a stable compoud.

## Prerequisites
Babel, shaep

### Convert xyz to mol2
$ ```babel -ixyz *.xyz -omol2 *.mol2```

### One molecule from original mol2 to restored from notation

```
from mol2_worker import atoms_and_bonds, xyz_names_bonds  
from mol2_chain import mol2_to_notation, dimensional_structure, write_mol2_file

atoms_info = atoms_and_bonds(name + '.mol2') #get information about atoms from mol2 - it's needed for write new mol2 later
notation_keeper = mol2_to_notation(xyz_names_bonds(name + '.mol2')) #get notation
dim_structure = dimensional_structure([notation_keeper[0], bonds_of_paired(notation_keeper[1])]) #return to dimentional xyz
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

Using shaep $ ```shaep -q original.mol2 after_restore.mol2 output_difference_file```
Or:
```
import sys, os, subprocess, shutil, tempfile
d = os.getcwd()
tmpdir = tempfile.mkdtemp()

compared_molecule2 = os.path.join(d, name + '.mol2')
compared_molecule = os.path.join(d, 'long.mol2')

sim_file = os.path.join(tmpdir, 'sim.txt')
subprocess.call([os.path.join(d, "shaep"), "-q", compared_molecule2, compared_molecule, sim_file],
                    stdout=subprocess.DEVNULL)
with open(sim_file, 'r') as f:
    f.readline()
    print(f.readline())
    #coeffs = f.readline().split()
    # coeffs = coeffs[2:6] #+ coeffs[8:]
    # print(coeffs)
    # print(np.linalg.norm(np.array(list(map(float, coeffs)))))
shutil.rmtree(tmpdir)
```



