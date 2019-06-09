import numpy as np
import os
import tempfile
from shutil import rmtree
from searcher_in_space import read_report
from mopac_worker import get_energy_of_xyz
import matplotlib.pyplot as plt

def read_approx_report(report_name, xyz):
    with open(report_name, 'r') as f:
        k = int(f.readline())
        f.readline()
        ens = []
        poss = []
        flag = True
        while flag:
            mol = []
            for _ in range(k):
                mol.append(f.readline().split())
            with open(xyz, 'w') as ff:
                ff.write(str(k)+'\n\n')
                for pos in mol:
                    ff.write('\t'.join([str(ix) for ix in pos])+'\n')
            ens.append(get_energy_of_xyz(xyz))
            flag = f.readline().replace(' ','').replace('\n', '').isdigit()
            poss.append(mol)
            if flag: f.readline()
    return ens

reports = '/home/anastasiia/PycharmProjects/chem_structure/mopac_html_example/js_example/ts_track.xyz'
# reports = ['equations_mnw_3a4_0'+str(i) for i in range(50)]
tmp = tempfile.mkdtemp()
xyz = os.path.join(tmp, 'coord.xyz')
mops = []
ens = read_approx_report(reports, xyz)

print(ens)
rmtree(tmp)

#
# for i in range(len(ch)):
#     if ch[i] == False:
#         with open(xyz, 'w') as ff:
#             ff.write(str(len(poss[i]))+'\n\n')
#             for pos in poss[i]:
#                 ff.write('\t'.join([str(ix) for ix in pos])+'\n')
#         appr_c.append(get_energy_of_xyz(xyz))
#         mops.append(ens[i])
# m = np.array(mops[:100])
# a = np.array(appr_c[:100])
# print(np.corrcoef(m,a))
# print(np.mean(m-a))
# plt.scatter(mops[:100:], appr_c[:100:])
# plt.xlabel('Mopac energy, eV')
# plt.ylabel('Approximate energy, eV')
# plt.show()

# print(max(ens))


#
# pred = []
# to_names = poss.pop(0)
# names = [i[0] for i in to_names]
# atoms = [np.array(list(map(float, i[1::]))) for i in to_names]
# ix = 0
# popper = []
# print(len(ens))
# for item, energy in zip(poss, ens):
#     atoms = [np.array(list(map(float, i[1::]))) for i in item]
#         if abs(pred[-1] - ens[ix]) > 2:
#             pred.pop(-1)
#             popper.append(ix)
#     else:
#         popper.append(ix)
#     ix += 1
#
# for i in popper[::-1]:
#     ens.pop(i)
# print(np.corrcoef(np.array(ens), np.array(pred)))
# print(len(ens))