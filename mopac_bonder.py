from subprocess import call
import os, glob, tempfile, shutil
from mol2_worker import *
from mopac_worker import *
# from mol2_chain import write_mol2_file

FNULL = open(os.devnull, 'w')
mopac_alias = '/opt/mopac/run_script.sh'
greece_numericals = {1: 'SINGLET', 2: 'DOUBLET',
                    3: 'TRIPLET', 4: 'QUARTET', 5: 'QUINTET',
                     6: 'SEXTET'}
calculations = {'Single Point': 'NOOPT',
                'Equilibrium Geometry': '',
                'Frequencies': 'FORCE'}


def generate_bond_finder_mop(file_xyz, file_to_bond):
    title = 'It is a bond finder file'
    calculate = 'NOOPT BONDS'
    theory = 'PM7'
    multiplicity, charge = 1, 0
    output = ''
    multStr = greece_numericals.get(multiplicity, 'SINGLET')

    # Calculation type:
    output += ' AUX LARGE CHARGE=%d %s %s %s\n' % \
              (charge, multStr, calculate, theory)
    output += '%s\n\n' % title
    positions = []
    with open(file_xyz, 'r') as f_r:
        n = int(f_r.readline())
        f_r.readline()
        for _ in range(n):
            positions.append(f_r.readline().split())

    for line in positions:
        output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(line[0], line[1], line[2], line[3])

    with open(file_to_bond, 'w') as f_w:
        f_w.write(output)
    return output


if __name__ == '__main__':
    tmpdir = tempfile.mkdtemp()
    example_xyz = './tmp/4b-Mn-OMe-ads-MeCHO_opt.xyz'
    with open(example_xyz, 'r') as f:
        n = int(f.readline())
    finder = os.path.join(tmpdir, '4b-Mn-OMe-ads-MeCHO_opt.mop')
    generate_bond_finder_mop(example_xyz, finder)
    call([mopac_alias, finder])
    with open(finder[:-3]+'out', 'r') as f:
        next(l for l in f if '(VALENCIES)   BOND ORDERS' in l)
        for i in range(n):
            # next(l for l in f if l!='')
            print(f.readline())
    shutil.rmtree(tmpdir)
