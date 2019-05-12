#!/bin/bash

import sys

file = sys.argv[1] #*.out after MOPAC run
lines = [line for line in open(file, 'r')]

n_c = lines.index('                             CARTESIAN COORDINATES\n')+2
count_c = lines[n_c::].index('\n')
lines = lines[n_c:n_c+count_c]
with open(sys.argv[2], 'w') as f:
    f.write(str(len(lines))+'\n')
    for i in lines:
        f.write('\n')
        f.write('\t'.join([syms for syms in i.split()[1::]]))