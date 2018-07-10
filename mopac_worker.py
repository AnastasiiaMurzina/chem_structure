from subprocess import call
import os
from mol2_worker import *
from mol2_chain import write_mol2_file

FNULL = open(os.devnull, 'w')
greece_numericals = {1: 'SINGLET', 2: 'DOUBLET',
                    3: 'TRIPLET', 4: 'QUARTET', 5: 'QUINTET',
                     6: 'SEXTET'}
calculations = {'Single Point': 'NOOPT',
                'Equilibrium Geometry': '',
                'Frequencies': 'FORCE'}

def generateInputFile(opts, mol_file_positions, mol_file_names):
    title = opts['Title']
    calculate = opts['Calculation Type']
    theory = opts['Theory']
    multiplicity = opts['Multiplicity']
    charge = opts['Charge']

    output = ''

    multStr = greece_numericals.get(multiplicity, 'SINGLET')

    # Calculation type:
    calcStr = calculations[calculate]

    # Charge, mult, calc type, theory:
    output += ' AUX LARGE CHARGE=%d %s %s %s\n' %\
        (charge, multStr, calcStr, theory)
    output += '%s\n\n' % title

    if calculate == 'Single Point':
        for num, pos in mol_file_positions.items():
            output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(mol_file_names[num].name, pos[0], pos[1], pos[2])
    else:
        for num, pos in mol_file_positions.items():
            output += ' {0}\t{1}\t1\t{2}\t1\t{3}\t1\n'.format(mol_file_names[num].name, pos[0], pos[1], pos[2])

    return output

def writeInputFile(opt, positions, names, file_name):
    with open(file_name, 'w') as f:
        f.write(generateInputFile(opt, positions, names))

def mopacOut_to_xyz(mopac_file, outXYZ_file):
    with open(mopac_file+'.arc', 'r') as f:
        for _ in range(8):
            f.readline()
        info = f.readline().split()
        count = int(info[-2])
        formula = (info[:-3])
        for _ in range(11):
            f.readline()
        energy = float(f.readline().split()[-2])
    with open(mopac_file + '.out', 'r') as f:
        next(l for l in f if 'CARTESIAN COORDINATES' in l)
        for _ in range(3):
            next(f)
        coords = [next(f) for _ in range(count)]
    with open(outXYZ_file, 'w') as f1:
        f1.write(str(count)+'\n')
        f1.write(' '.join(formula)+'\n')
        for i in coords:
            line = i.split()
            f1.write('{}\t{}\t{}\t{}\n'.format(line[1], line[2], line[3], line[4]))
    return energy


if __name__ == '__main__':
    ''' Compare optimizations with Cartesian and Z-matrix'''
    options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Equilibrium Geometry',
                              'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM6'}
    # writeInputFile(options, )

    # file1 = 'Aniline_cart_av'
    # file2 = 'Aniline_zmax_av'
    # energy1 = mopacOut_to_xyz(file1, file1 + '_opt.xyz')
    # energy2 = mopacOut_to_xyz(file2, file2 + '_opt.xyz')
    # print(energy1, energy2)
    # call(['babel', '-ixyz', file1+'_opt.xyz', '-omol', file1 + '_opt.mol'])
    # call(['babel', '-ixyz', file2+'_opt.xyz', '-omol', file2 + '_opt.mol'])
    # call(['./shaep', '-q', file1 + '_opt.mol', file2 + '_opt.mol', 'output'])



    # file_name = 'My_Aniline.mol2'
    # positions = xyz_names(file_name)
    # bonds, names = xyz_names_bonds(file_name)
    # atom_info = atoms_and_bonds(file_name)
    # mopac_file_name = 'My_Aniline'

    # writeInputFile(options, positions, names, mopac_file_name+'.mop')
    # call(['/opt/mopac/run_script.sh', mopac_file_name + '.mop'],
    #             stdout=FNULL)
    # energy = mopacOut_to_xyz(mopac_file_name, 'My_Aniline_opt.xyz')
    # print(energy)
    pass



