import sys, os, subprocess, shutil, tempfile
import numpy as np
from itertools import product
from mol2_chain import atoms_and_bonds, mol2_to_notation, xyz_names_bonds, bonds_of_paired, dimensional_structure, write_mol2_file, xyz_names
from mopac_worker import writeInputFile, mopacOut_to_xyz, writeInputFileFromXYZ, mopacOut_to_xyz_with_energy

name = 'Mn-deprot-Mn-bare'

FNULL = open(os.devnull, 'w')
# mopac_alias = 'mopac'
mopac_alias = '/opt/mopac/run_script.sh'
tmpdir = tempfile.mkdtemp()
initial_dir ='../mols_dir'
here_dir = os.getcwd()
options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
                   'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}

initial_xyz_many = [f for f in os.listdir(initial_dir) if os.path.isfile(os.path.join(initial_dir, f))]
for one in initial_xyz_many:
    clear_name = one.split('.')[0]
    first_mop = os.path.join(tmpdir, clear_name+'.mop')
    first_mop_name = os.path.join(tmpdir, clear_name)
    names, positions = [], []
    writeInputFileFromXYZ(options, os.path.join(initial_dir, one), first_mop)
    subprocess.call([mopac_alias, first_mop], stdout=FNULL)
    mopacOut_to_xyz_with_energy(first_mop_name, clear_name+'_opt.xyz')
shutil.rmtree(tmpdir)
