import sys, os, subprocess, shutil, tempfile
import numpy as np
from berny2.berny import Berny, geomlib, optimize
from berny2.berny.solvers import MopacSolver
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file, xyz_names
from mopac_worker import writeInputFile, mopacOut_to_xyz
import time
names = ['Aniline', 'Caffein', 'Cocaine', 'Ethanol', 'Ethene']#, 'Heroin']
# names = ['Aniline', 'Caffein', 'Caffein_1', 'Cocaine', 'Ethanol', 'Ethene', 'Heroin', 'Heroin_1', 'Heroin_2',
#          'Naphthalene', 'Phenol', 'pyridine', 'Toluene']
sx = 0.05
# sx = (1+np.cos(0.4*np.pi))**0.5
times_and_currency = {}
FNULL = open(os.devnull, 'w')
# mopac_alias = 'mopac'
mopac_alias = '/opt/mopac/run_script.sh'
tmpdir = tempfile.mkdtemp() #make temporaly directory
options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
                       'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}

kcal, ev, angstrom = 627.503, 27.2107, 0.52917721092

for name in names:
    start_time = time.time()
    original_xyz = os.path.join(os.curdir, 'mols_dir', name + '.xyz') #get structure from ./mols_dir/
    mol_initial = geomlib.readfile(original_xyz)
    solver = MopacSolver()
    final = optimize(solver, mol_initial, steprms=0.01, stepmax=sx)
    original_opt_xyz = os.path.join(tmpdir, name+'_opt.xyz')
    with open(original_opt_xyz, 'w') as f:
        final.dump(f, 'xyz')

    original_mol2 = os.path.join(tmpdir, name+'.mol2')
    subprocess.call(['babel', '-ixyz', original_opt_xyz, '-omol2', original_mol2], stdout=FNULL)
    positions = xyz_names(original_mol2)
    bonds, names = xyz_names_bonds(original_mol2)
    atom_info = atoms_and_bonds(original_mol2)

    original_opt_mol2 = os.path.join(tmpdir, name +'_opt.mol2')
    original_opt_mol = os.path.join(tmpdir, name + '_opt.mol')
    subprocess.call(['babel', '-ixyz', original_opt_xyz, '-omol2', original_opt_mol2])
    subprocess.call(['babel', '-ixyz', original_opt_xyz, '-omol', original_opt_mol])

    atoms_info = atoms_and_bonds(original_opt_mol2) #get information about atoms from mol2 - it's needed for write new mol2 later
    notation_keeper = mol2_to_notation(xyz_names_bonds(original_opt_mol2), method='ten') #get notation
    paired = bonds_of_paired(notation_keeper[1])
    dim_structure = dimensional_structure([notation_keeper[0], paired],  method='ten') #return to dimentional xyz
    reconstructed_mol2 = os.path.join(tmpdir, 'my_'+name+'.mol2')
    # reconstructed_mol = os.path.join(tmpdir, 'my_'+name+'.mol')
    write_mol2_file(reconstructed_mol2, atoms_info, dim_structure, bonds=paired) #write restored dimentional structure
    reconstructed_xyz = os.path.join(tmpdir, 'my_'+name+'.xyz')
    subprocess.call(['babel', '-imol2', reconstructed_mol2, '-oxyz', reconstructed_xyz], stdout=FNULL)
    reconstructed_opt_xyz = os.path.join(tmpdir, 'my_' + name + '_opt.xyz')
    reconstructed_opt_mol = os.path.join(tmpdir, 'my_' + name + '_opt.mol')
    reconstructed_opt_mol2 = os.path.join(tmpdir, 'my_' + name + '_opt.mol2')

    mol_reconstructed = geomlib.readfile(reconstructed_xyz)
    solver_rec = MopacSolver()
    final_rec = optimize(solver_rec, mol_reconstructed, steprms=0.01, stepmax=sx)
    with open(reconstructed_opt_xyz, 'w') as f:
        final_rec.dump(f, 'xyz')

    subprocess.call(['babel', '-ixyz', reconstructed_opt_xyz, '-omol', reconstructed_opt_mol])
    subprocess.call(['babel', '-ixyz', reconstructed_opt_xyz, '-omol2', reconstructed_opt_mol2])
    name_check = os.path.join(tmpdir, name)
    subprocess.call(['./shaep', '-q', original_opt_mol, reconstructed_opt_mol, name_check])#, ' --transformDistance=1.9'], stdout=FNULL)
    with open(name_check, 'r') as f:
        for line in f:
            pass
        print(line)
    working_time = time.time() - start_time

    mopac_input = '{} 1SCF GRADIENTS\n\n\n'.format(options['Theory']) + '\n'.join('{} {} 1 {} 1 {} 1'.format(el, *coord) for el, coord in geomlib.readfile(original_opt_xyz))
    input_file = os.path.join(tmpdir, 'job.mop')
    with open(input_file, 'w') as f:
        f.write(mopac_input)
    subprocess.check_output(['/opt/mopac/run_script.sh', input_file])
    with open(os.path.join(tmpdir, 'job.out')) as f:
        energy = float(next(l for l in f if 'TOTAL ENERGY' in l).split()[3])/ev

    mopac_input = '{} 1SCF GRADIENTS\n\n\n'.format(options['Theory']) + '\n'.join('{} {} 1 {} 1 {} 1'.format(el, *coord) for el, coord in geomlib.readfile(reconstructed_opt_xyz)
                )
    input_file2 = os.path.join(tmpdir, 'job2.mop')
    with open(input_file2, 'w') as f:
        f.write(mopac_input)
    subprocess.check_output(['/opt/mopac/run_script.sh', input_file])
    with open(os.path.join(tmpdir, 'job.out')) as f:
        energy2 = float(next(l for l in f if 'TOTAL ENERGY' in l).split()[3])/ev

    times_and_currency.update({name: [working_time, energy, energy2, line]})

    # print(final.current, final_rec.current)
shutil.rmtree(tmpdir)
print(times_and_currency)
