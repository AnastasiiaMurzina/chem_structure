#!/bin/bash

import sys

file = sys.argv[1]
lines = [line for line in open(file, 'r')]

with open(file, 'w') as f:
    for i in lines:
        f.write('\n')
        f.write('\t'.join([syms for syms in i.split()[1::]]))