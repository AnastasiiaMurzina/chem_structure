import sys, os, subprocess, shutil, tempfile
import numpy as np
from itertools import product
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file, xyz_names
from mopac_worker import writeInputFile, mopacOut_to_xyz

name = 'Caffein'

FNULL = open(os.devnull, 'w')
mopac_alias = '/opt/mopac/run_script.sh'
tmpdir = tempfile.mkdtemp()
initial_xyz = os.path.join(os.curdir, 'mols_dir', name+'.xyz')

options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
                   'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM6'}
original_mol2 = os.path.join(tmpdir, name+'.mol2')
subprocess.call(['babel', '-ixyz', initial_xyz, '-omol2', original_mol2], stdout=FNULL)
positions = xyz_names(original_mol2)
bonds, names = xyz_names_bonds(original_mol2)
atom_info = atoms_and_bonds(original_mol2)
first_mop = os.path.join(tmpdir, name+'.mop')
first_mop_name = os.path.join(tmpdir, name)
writeInputFile(options, positions, names, first_mop)
subprocess.call([mopac_alias, first_mop], stdout=FNULL)

original_opt_xyz = os.path.join(tmpdir, name+'_opt.xyz')

energy1 = mopacOut_to_xyz(first_mop_name, original_opt_xyz)

original_opt_mol2 = os.path.join(tmpdir, name+'_opt.mol2')
original_opt_mol = os.path.join(tmpdir, name+'_opt.mol')
subprocess.call(['babel', '-ixyz', original_opt_xyz, '-omol2', original_opt_mol2])
subprocess.call(['babel', '-ixyz', original_opt_xyz, '-omol', original_opt_mol])

atoms_info = atoms_and_bonds(original_opt_mol2) #get information about atoms from mol2 - it's needed for write new mol2 later
notation_keeper = mol2_to_notation(xyz_names_bonds(original_opt_mol2), method='ten') #get notation
paired = bonds_of_paired(notation_keeper[1])
dim_structure = dimensional_structure([notation_keeper[0], paired], method='ten') #return to dimentional xyz
reconstructed_mol2 = os.path.join(tmpdir, 'my_'+name+'.mol2')
reconstructed_mol = os.path.join(tmpdir, 'my_'+name+'.mol')
write_mol2_file(reconstructed_mol2, atoms_info, dim_structure, bonds=paired) #write restored dimentional structure

reconstructed_mop = os.path.join(tmpdir, 'my_' + name + '.mop')
reconstructed_mop_name = os.path.join(tmpdir, name)
writeInputFile(options, positions, names, reconstructed_mop)
subprocess.call([mopac_alias, reconstructed_mop], stdout=FNULL)
reconstructed_opt_xyz = os.path.join(tmpdir, 'my_' + name + '_opt.xyz')
reconstructed_opt_mol = os.path.join(tmpdir, 'my_' + name + '_opt.mol')
reconstructed_opt_mol2 = os.path.join(tmpdir, 'my_' + name + '_opt.mol2')
energy2 = mopacOut_to_xyz(reconstructed_mop_name, reconstructed_opt_xyz)
subprocess.call(['babel', '-ixyz', reconstructed_opt_xyz, '-omol', reconstructed_opt_mol])
subprocess.call(['babel', '-ixyz', reconstructed_opt_xyz, '-omol2', reconstructed_opt_mol2])
print(energy1, energy2)
subprocess.call(['./shaep', '-q', original_opt_mol, reconstructed_opt_mol, os.path.join(tmpdir,name)], stdout=FNULL)
with open(os.path.join(tmpdir, name), 'r') as f:
    for line in f:
        pass
    print(line)
# subprocess.call(['babel', '-imol2'])

# m2_mol = os.path.join(tmpdir, 'my_'+name+'.mol')
# subprocess.call(["babel", "-imol2", m2_mol2, '-omol', m2_mol],
#                 stdout=FNULL)
# m2_mol_opt = os.path.join(tmpdir, 'my_'+name+'_opt.mol')

shutil.rmtree(tmpdir)