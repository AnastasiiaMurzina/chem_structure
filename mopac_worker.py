from subprocess import call
import os, glob, tempfile, shutil

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
    if opts.get('TS'): # Locate-ts, saddle???
        output += ' TS\n'
    else:
        output += ' AUX LARGE CHARGE=%d %s %s %s\n' %\
        (charge, multStr, calcStr, theory)
    output += '%s\n\n' % title

    if calculate == 'Single Point':
        for num, pos in mol_file_positions.items():
            output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(mol_file_names[num].name, *pos[0:3:])
    else:
        for num, pos in mol_file_positions.items():
            output += ' {0}\t{1}\t1\t{2}\t1\t{3}\t1\n'.format(mol_file_names[num].name, *pos[0:3:])

    return output


def writeInputFileFromXYZ(opts, xyz_file, mop_file):
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
    if opts.get('TS'):
        output += ' TS\n'
    else:
        output += ' AUX LARGE CHARGE=%d %s %s %s\n' %\
        (charge, multStr, calcStr, theory)
    output += '%s\n\n' % title
    with open(xyz_file, 'r') as f:
        n = int(f.readline())
        f.readline()
        lines = [f.readline().split() for _ in range(n)]

    if calculate == 'Single Point':
        for line in lines:
            output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(line[0], line[1], line[2], line[3])
    else:
        for line in lines:
            output += ' {0}\t{1}\t1\t{2}\t1\t{3}\t1\n'.format(line[0], line[1], line[2], line[3])
    with open(mop_file, 'w') as f:
        f.write(output)
    return output


def writeInputFileFromMOL2(opts, mol2_file, mop_file):
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
    if opts.get('TS'):
        output += ' AUX LARGE CHARGE=%d %s %s %s %s\n' %\
        (charge, multStr, calcStr, theory, opts['TS'])
    else:
        output += ' AUX LARGE CHARGE=%d %s %s %s\n' %\
        (charge, multStr, calcStr, theory)
    output += '%s\n\n' % title
    with open(mol2_file, 'r') as f:
        next(l for l in f if '@<TRIPOS>ATOM' in l)
        s = f.readline()
        coords = []
        while '@<TRIPOS>BOND' not in s:
            coords.append(s.split()[1:5:])
            s = f.readline()
    if calculate == 'Single Point':
        for line in coords:
            output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(line[0], line[1], line[2], line[3])
    else:
        for line in coords:
            output += ' {0}\t{1}\t1\t{2}\t1\t{3}\t1\n'.format(line[0], line[1], line[2], line[3])
    with open(mop_file, 'w') as f:
        f.write(output)
    return output


def writeInputFile(opt, positions, names, file_name):
    with open(file_name, 'w') as f:
        f.write(generateInputFile(opt, positions, names))


def mopacOut_to_xyz(mopac_file, outXYZ_file):
    '''
    :param mopac_file: file name without extension performed mop-process
    :param outXYZ_file: xyz to write file
    :return: None
    '''
    with open(mopac_file+'.arc', 'r') as farc:
        for _ in range(8):
            farc.readline()
        info = farc.readline().split()
        count = int(info[-2])
        formula = (info[:-3])
    with open(mopac_file + '.out', 'r') as f:
        next(l for l in f if 'CARTESIAN COORDINATES' in l)
        next(l for l in f if 'CARTESIAN COORDINATES' in l)
        next(f)
        coords = [next(f) for _ in range(count)]
    with open(outXYZ_file, 'w') as f1:
        f1.write(str(count)+'\n')
        f1.write(' '.join(formula)+'\n')
        for i in coords:
            line = i.split()
            f1.write('{}\t{}\t{}\t{}\n'.format(line[1], line[2], line[3], line[4]))


def get_energy(mopac_file):
    with open(mopac_file, 'r') as f:
        try:
            line = next(l for l in f if 'TOTAL ENERGY' in l)
        except:
            return None
        energy = float(line.split()[-2])
    return energy


def get_heat(mopac_file):
    with open(mopac_file, 'r') as f:
        try:
            line = next(l for l in f if 'FINAL HEAT OF FORMATION' in l)
        except:
            return None
        return float(line.split()[5])


def get_energy_of_xyz(xyz_file, tmpdir='', devnull=True):
    del_flag = False
    if tmpdir == '':
        tmpdir = tempfile.mkdtemp()
        del_flag = True
        file_name = os.path.basename(os.path.normpath(xyz_file))
    xyz_to_mop = os.path.join(tmpdir, file_name[:-4:] + '.mop')
    header = ' AUX LARGE CHARGE=0 NOOPT PM7\nTitle\n'
    with open(xyz_file, 'r') as f:
        with open(os.path.join(tmpdir, xyz_to_mop), 'w') as f_w:
            f_w.write(header)
            n = int(f.readline())
            f.readline()
            for _ in range(n):
                line = f.readline().split()
                f_w.write('{}\t{}\t0\t{}\t0\t{}\t0\n'.format(*line))
    call(['mopac', os.path.join(tmpdir, xyz_to_mop)], stdout=FNULL, stderr=FNULL)
    a = get_energy(os.path.join(tmpdir, xyz_to_mop[:-4])+'.out')
    if del_flag: shutil.rmtree(tmpdir)
    return a


def get_heat_of_xyz(xyz_file, tmpdir=''):
    del_flag = False
    if tmpdir == '':
        tmpdir = tempfile.mkdtemp()
        del_flag = True
        file_name = os.path.basename(os.path.normpath(xyz_file))
    xyz_to_mop = os.path.join(tmpdir, file_name  + '.mop')
    header = ' AUX LARGE CHARGE=0 NOOPT PM7\nTitle\n'
    with open(xyz_file, 'r') as f:
        with open(os.path.join(tmpdir, xyz_to_mop), 'w') as f_w:
            f_w.write(header)
            n = int(f.readline())
            f.readline()
            for _ in range(n):
                line = f.readline().split()
                f_w.write('{}\t{}\t0\t{}\t0\t{}\t0\n'.format(*line))
    call(['run_mopac', os.path.join(tmpdir, xyz_to_mop)])
    a = get_heat(os.path.join(tmpdir, file_name )+'.out')
    if del_flag: shutil.rmtree(tmpdir)
    return a


def get_energy_of_mol2(mol2_file):
    tmpdir = tempfile.mkdtemp()
    file_name = os.path.basename(os.path.normpath(mol2_file))
    xyz_to_mop = os.path.join(tmpdir, file_name[:-5:] + '.xyz')
    call(['babel', '-imol2', mol2_file, '-oxyz', xyz_to_mop])
    return get_energy_of_xyz(xyz_to_mop, tmpdir=tmpdir)


def get_xyz_files(directory):
    return glob.glob(directory+'/*.xyz')


def get_mop_files(directory):
    return glob.glob(directory+'/*.mop')



if __name__ == '__main__':
    ''' Compare optimizations with Cartesian and Z-matrix'''
    options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Single Point',
                              'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}

    def run_mopac():
        mopac_alias = 'mopac'
        tmpdir = tempfile.mkdtemp()
        g_dir = os.path.join(os.getcwd(), 'q_compare')
        for i in [x[0] for x in os.walk(g_dir)]:
            files = glob.glob(i+'/*.mol2')
            for j in files:
                name = os.path.join(tmpdir, j[:-5])
                writeInputFileFromMOL2(options, j, name+'.mop')
                call([mopac_alias, name+'.mop'])
                mopacOut_to_xyz(name, os.path.join(os.getcwd(), 'q_compare', j[:-5:] + '_opt.xyz'))
        shutil.rmtree(tmpdir)

    for i in ['2-Mn-OtBu_opt.mol2', '3-MnH2_opt.mol2', '4c-Mn-OMe_opt.mol2']:
        name = os.path.join('./tmp', i)
        print(i, get_energy_of_mol2(name))


    def get_opt():
        g_dir = os.path.join(os.getcwd(), 'mols_dir', 'Molecules')
        for i in [x[0] for x in os.walk(g_dir)][1::]:
            cur_d = os.path.join(g_dir, i.split('/')[-1])
            files = glob.glob(i+'/*.arc')
            for j in files:
                mopacOut_to_xyz(j[:-4:], os.path.join(cur_d, j[:-4:]+'_opt.xyz'))

    pass



