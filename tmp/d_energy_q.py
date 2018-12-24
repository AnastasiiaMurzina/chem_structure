import os, subprocess, tempfile, shutil
from mol2_chain_q import mol2_to_notation, dimensional_structure, write_mol2_file, bonds_of_paired, atoms_and_bonds, xyz_names_bonds
from mopac_worker import writeInputFileFromXYZ,  mopacOut_to_xyz_with_energy

FNULL = open(os.devnull, 'w')
mopac_alias = '/opt/mopac/run_script.sh'
options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Single Point',
                              'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}

if __name__ == '__main__':
    tmpdir = tempfile.mkdtemp()
    files = [f for f in os.listdir(os.getcwd()) if os.path.isfile(os.path.join(os.getcwd(), f))]
    for file in files:
        name = file.split('.')[-2]
        mol2_original = os.path.join(tmpdir, name+'.mol2')
        subprocess.call(['babel', '-ixyz', file, '-omol2', mol2_original])

        atoms_info = atoms_and_bonds(mol2_original)
        ln = mol2_to_notation(xyz_names_bonds(mol2_original))
        paired = bonds_of_paired(ln[1])
        dim_structure = dimensional_structure([ln[0], paired])

        br_mol2 = os.path.join(tmpdir, 'My_q_'+name+'.mol2')
        br_xyz = os.path.join(tmpdir, 'My_q_'+name+'.xyz')
        write_mol2_file(br_mol2, atoms_info, dim_structure, bonds=paired)
        subprocess.call(['babel', '-imol2', br_mol2, '-oxyz', br_xyz])
        mop_file = os.path.join(tmpdir, name+'.mop')
        mop_file_cl = os.path.join(tmpdir, name)

        writeInputFileFromXYZ(options, br_xyz, mop_file)
        subprocess.call([mopac_alias, mop_file], stdout=FNULL)
        en1 = mopacOut_to_xyz_with_energy(mop_file_cl, os.path.join(tmpdir, name + '_opt.xyz'))
        with open(file, 'r') as f:
            f.readline()
            en0 = float(f.readline())
        print(name, en0, en1, abs(en0-en1), abs((en0-en1)/en0))
    shutil.rmtree(tmpdir)