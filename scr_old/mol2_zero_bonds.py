# import tempfile
import subprocess

if __name__ == '__main__':
    # tmpf = tempfile.TemporaryFile()
    tmpf = open('TS-2c-2d.mop', 'w')
    tmpf.write('BONDS\nTitle\n\n')
    with open('TS-2c-2d_last.xyz', 'r') as f:
        n = int(f.readline())
        f.readline()
        for _ in range(n):
            line = f.readline().split()
            tmpf.write("{}\t{}\t0\t{}\t0\t{}\t0\n".format(*line))
    tmpf.close()
    subprocess.call(['/opt/mopac/run_script.sh', 'TS-2c-2d.mop'])

    pass




