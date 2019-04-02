from os import path
from tempfile import mkdtemp
from shutil import rmtree
from mol2_worker import xyz_to_array
from subprocess import call, Popen, PIPE


def get_dG(file_xyz, tmp=''):
    flag_deleter = False
    if tmp == '':
        flag_deleter = True
        tmp = mkdtemp()

    just_name = path.basename(file_xyz)[:-4:]
    product = xyz_to_array(file_xyz)
    mop = path.join(tmp, just_name + '.mop')
    out = path.join(tmp, just_name + '.out')
    with open(mop, 'w') as f_w:
        f_w.write('LET CHARGE=0 SINGLET PM6-D3H4X NOOPT THERMO\nTi_C_complex_UFF\n\n')
        for i in product:
            f_w.write('{}\t{}\t{}\t{}\n'.format(i[0], str(i[1]), str(i[2]), str(i[3])))
    call(['/opt/mopac2/run_mopac', mop])
    # proc = Popen(['sudo', '-S', '/opt/mopac/run_script.sh', mop], stdout=PIPE, stdin=PIPE, stderr=PIPE, universal_newlines=True)
    # proc.stdin.write("\n")
    # out, err = proc.communicate(input="\n")

    with open(out, 'r') as f:
        next(l for l in f if 'HEAT OF FORMATION' in l)
        next(l for l in f if 'HEAT OF FORMATION' in l)
        f.readline()
        zpe_line = f.readline().split()
        zpe = float(zpe_line[3])
        next(l for l in f if 'CALCULATED THERMODYNAMIC PROPERTIES' in l)
        line = f.readline().split()
        while line[0] != '370.00':
            line = f.readline().split()
            while line == []:
                line = f.readline().split()
        for _ in range(3): f.readline()
        line = f.readline().split()
        entalpy, entropy = line[2], line[4]
        entalpy, entropy = float(entalpy), float(entropy)
        dG = entalpy - 370 * entropy

    if flag_deleter: rmtree(tmp)
    return entropy, entalpy, dG, zpe


if __name__ == '__main__':
    pass