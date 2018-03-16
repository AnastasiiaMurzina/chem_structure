from chain_with_rotate import read_xyz_file, coordinates_to_notation
from itertools import combinations
import numpy as np
# coords, names = read_xyz_file('pyridine.xyz')
# print(len(names))
# for i, j in combinations(range(1,len(names)+1), 2):
#     print(i,j)
#     print(np.linalg.norm(coords[i]-coords[j]))
print(coordinates_to_notation(read_xyz_file('Caffein.xyz'), valid_length=0.7))