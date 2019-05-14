from subprocess import call
import os, glob, tempfile, shutil
# from mol2_worker import *
# from mol2_chain import write_mol2_file

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

    # with open(mopac_file+'.out', 'r') as f:
    #     for _ in range(8):
    #         f.readline()
    #     info = f.readline().split()

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
            next(l for l in f if 'FINAL HEAT OF FORMATION' in l)
        except:
            return None
        next(f)
        next(f)
        energy = float(next(f).split()[-2])
    return energy

def get_heat(mopac_file):
    with open(mopac_file, 'r') as f:
        try:
            line = next(l for l in f if 'FINAL HEAT OF FORMATION' in l)

        except:
            return None
        return float(line.split()[5])


def get_energy_of_xyz(xyz_file, tmpdir=''):
    del_flag = False
    if tmpdir == '':
        tmpdir = tempfile.mkdtemp()
        del_flag = True
    name = os.path.basename(os.path.normpath(xyz_file))
    xyz_to_mop = os.path.join(tmpdir, name[:-4:] + '.mop')
    header = ' AUX LARGE CHARGE=0 SINGLET NOOPT PM7\nTitle\n'
    with open(xyz_file, 'r') as f:
        with open(os.path.join(tmpdir, xyz_to_mop), 'w') as f_w:
            f_w.write(header)
            n = int(f.readline())
            f.readline()
            for _ in range(n):
                line = f.readline().split()
                f_w.write('{}\t{}\t0\t{}\t0\t{}\t0\n'.format(*line))
    call(['/opt/mopac2/run_mopac', os.path.join(tmpdir, xyz_to_mop)])
    a = get_energy(os.path.join(tmpdir, xyz_to_mop)[:-4]+'.out')
    if del_flag: shutil.rmtree(tmpdir)
    return a


def get_heat_of_xyz(xyz_file, tmpdir=''):
    del_flag = False
    if tmpdir == '':
        tmpdir = tempfile.mkdtemp()
        del_flag = True
    name = os.path.basename(os.path.normpath(xyz_file))
    xyz_to_mop = os.path.join(tmpdir, name + '.mop')
    header = ' AUX LARGE CHARGE=0 SINGLET NOOPT PM7\nTitle\n'
    with open(xyz_file, 'r') as f:
        with open(os.path.join(tmpdir, xyz_to_mop), 'w') as f_w:
            f_w.write(header)
            n = int(f.readline())
            f.readline()
            for _ in range(n):
                line = f.readline().split()
                f_w.write('{}\t{}\t0\t{}\t0\t{}\t0\n'.format(*line))
    call(['/opt/mopac2/run_mopac', os.path.join(tmpdir, xyz_to_mop)])
    a = get_heat(os.path.join(tmpdir, name)+'.out')
    if del_flag: shutil.rmtree(tmpdir)
    return a


def get_energy_of_mol2(mol2_file):
    tmpdir = tempfile.mkdtemp()
    name = os.path.basename(os.path.normpath(mol2_file))
    xyz_to_mop = os.path.join(tmpdir, name[:-5:] + '.xyz')
    call(['babel', '-imol2', mol2_file, '-oxyz', xyz_to_mop])
    return get_energy_of_xyz(xyz_to_mop, tmpdir=tmpdir)


def get_xyz_files(dir):
    return glob.glob(dir+'/*.xyz')


def get_mop_files(dir):
    return glob.glob(dir+'/*.mop')



if __name__ == '__main__':
    ''' Compare optimizations with Cartesian and Z-matrix'''
    options = {'Title': 'Smth info about optimization', 'Calculation Type': 'Single Point',
                              'Charge': 0, 'Multiplicity': 1, 'Theory': 'PM7'}


    import subprocess
    # print(get_energy_of_mol2('./q_compare/2_no_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/2_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/3_no_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/3_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/4_no_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/4_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/5_no_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./q_compare/5_relax/My_2-Mn-OtBu_opt_q.mol2'))
    # print(get_energy_of_mol2('./tmp/My_2-Mn-OtBu_opt_q.mol2'))
    def run_mopac():
        mopac_alias = '/opt/mopac/run_script.sh'
        tmpdir = tempfile.mkdtemp()
        g_dir = os.path.join(os.getcwd(), 'q_compare')
        for i in [x[0] for x in os.walk(g_dir)]:
            files = glob.glob(i+'/*.mol2')
            for j in files:
                name = os.path.join(tmpdir, j[:-5])
                writeInputFileFromMOL2(options, j, name+'.mop')
                subprocess.call([mopac_alias, name+'.mop'])
                mopacOut_to_xyz(name, os.path.join(os.getcwd(), 'q_compare', j[:-5:] + '_opt.xyz'))
        shutil.rmtree(tmpdir)
            # files = get_mop_files(i)
            # for j in files:
            #     subprocess.call([mopac_alias, j])
    for i in ['2-Mn-OtBu_opt.mol2', '3-MnH2_opt.mol2', '4c-Mn-OMe_opt.mol2']:
        name = os.path.join('./tmp', i)
        print(i, get_energy_of_mol2(name))
    # g_dir = os.path.join(os.getcwd(), 'q_compare')
    # name = '2-Mn-OtBu_opt_q'
    # for i in [x[0] for x in os.walk(g_dir)]:
    #     files = glob.glob(i+'/*_opt.xyz')
    #     for j in files:
    #         print(j, get_energy_of_xyz(j))

    # run_mopac()

    def get_opt():
        g_dir = os.path.join(os.getcwd(), 'mols_dir', 'Molecules')
        # g_dir2 = os.path.join(os.getcwd(), 'mols_dir', 'Molecules_Mopt')
        for i in [x[0] for x in os.walk(g_dir)][1::]:
            cur_d = os.path.join(g_dir, i.split('/')[-1])
            # os.mkdir(cur_d)
            files = glob.glob(i+'/*.arc')
            for j in files:
                mopacOut_to_xyz(j[:-4:], os.path.join(cur_d, j[:-4:]+'_opt.xyz'))

    # get_opt()
    # g_dir = os.path.join(os.getcwd(), 'mols_dir', 'Molecules')
    # print([x[0] for x in os.walk(g_dir)][30])
    # k = 0
    # for i in [x[0] for x in os.walk(g_dir)]:
    #     k+=1
    #     print(k, i)

            # with open(j, 'r') as f:
            #     n = int(f.readline())
            #     f.readline()
            #     lines = []
            #     for _ in range(n):
            #         lines.append(f.readline().split())
            # title = options['Title']
            # calculate = options['Calculation Type']
            # theory = options['Theory']
            # multiplicity = options['Multiplicity']
            # charge = options['Charge']
            # output = ''
            # multStr = greece_numericals.get(multiplicity, 'SINGLET')
            # calcStr = calculations[calculate]
            # if options.get('TS'):
            #     output += ' AUX LARGE CHARGE=%d %s %s %s %s\n' % \
            #               (charge, multStr, calcStr, theory, options['TS'])
            # else:
            #     output += ' AUX LARGE CHARGE=%d %s %s %s\n' % \
            #               (charge, multStr, calcStr, theory)
            # output += '%s\n\n' % title
            #
            # if calculate == 'Single Point':
            #     for i in lines:
            #         output += '  {0}\t{1}\t0\t{2}\t0\t{3}\t0\n'.format(i[0], i[1], i[2], i[3])
            # print(os.path.join(tmpdir, j[:-3:]+'mop'))
            # with open(os.path.join(tmpdir, j[:-3:]+'mop'), 'w') as f_w:
            #
            #     f_w.writelines(output)


    # writeInputFile(options, )
    # mopacOut_to_xyz('/home/anastasiia/PycharmProjects/chem_structure/tmp/TS-2b-2c_step1', '/home/anastasiia/PycharmProjects/chem_structure/tmp/TS-2b-2c_step1_from_out.xyz')

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



