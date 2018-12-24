n = 0
atoms = 55
name = 'TS-2b-2c OtBu-step1'
with open(name + '.log', 'r') as f:
    n = f.read().count('Coordinates (')
# print(n)
with open(name+'.log', 'r') as f:
    for i in range(n):
        name_c = name + str(i+1) + '.xyz'
        with open(name_c, 'w') as f_w:
            (next(l for l in f if 'Coordinates (' in l))
            next(f)
            next(f)
            f_w.write(str(atoms))
            f_w.write('\n\n')
            for _ in range(atoms):
                line = next(f).split()
                f_w.write('\t'.join([line[1], line[3], line[4], line[5]]))
                f_w.write('\n')

# def coords_from_gaussian_to_file(gaussian_file, xyz_file='', num=0):
#     with open(gaussian_file, 'r') as f_r:
#         print(f_r.count('Cartesian ('))