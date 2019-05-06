from os import path
import mol2_worker
from mopac_worker import get_energy_of_xyz, get_heat_of_xyz
from thermo_characteristics import get_dG
import tempfile, shutil
from thermo_characteristics import get_dG
import numpy as np
import matplotlib.pyplot as plt
import rmsd

def heats_from_path(file_xyz):
    states = []
    n = 24
    x, y = list(range(199)), []
    tmp = tempfile.mkdtemp()
    tmp_file = path.join(tmp, 'current.xyz')
    with open(file_xyz, 'r') as f:
        for _ in range(199):
            f.readline()
            f.readline()
            states += [f.readline().split() for _ in range(n)]
        for indx in x:
            current_state = states[indx*24:(indx+1)*24]
            with open(tmp_file, 'w') as f:
                f.write(str(n) + '\n\n')
                for inx in range(n):
                    f.write('{}\t{}\t{}\t{}\n'.format(*current_state[inx]))
            y.append(get_heat_of_xyz(tmp_file, tmpdir=tmp))
        shutil.rmtree(tmp)
    return x, y

# plt.plot(x, y)
# plt.title('Total heat from xyz (probably irc path)')
# plt.ylabel('KCAL/MOL')
# plt.show()