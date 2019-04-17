from math import pi, cos, sin
import rmsd
import numpy as np

file1 = './opted_mols/2-Mn-OtBu_opted.xyz'
with open(file1, 'r') as f1:
    n = int(f1.readline())
    f1.readline()
    coord = []
    for _ in range(n):
        coord.append(list(map(float, f1.readline().split()[1::])))

coord = np.array(coord)
coord2 = np.array([c+np.array([1, 1, 4]) for c in coord])

coord -= rmsd.centroid(coord)
coord2 -= rmsd.centroid(coord2)
print(rmsd.rmsd(coord, coord2))

phi, psi = pi/10, pi/4
coord3 = np.array([np.array([c[0]*(1), c[1]*cos(psi)+c[2]*sin(psi),
                             -c[1]*sin(psi)+c[2]*cos(psi)]) for c in coord])

coord3 -= rmsd.centroid(coord3)
rotate = rmsd.kabsch(coord3, coord)

coord3 = np.dot(coord3, rotate)

print(rmsd.rmsd(coord, coord3))
