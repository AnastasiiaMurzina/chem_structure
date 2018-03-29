import numpy as np

def read_mol2(file_name):
    with open(file_name, 'r') as f1:
        all_f = f1.read()
    ix0 = all_f.index('@<TRIPOS>MOLECULE') # info = [ix0:ix1]
    ix1 = all_f.index('@<TRIPOS>ATOM')     # xyz_part = [ix1:ix2]
    ix2 = all_f.index('@<TRIPOS>BOND')     # bond = [ix2::]
    return all_f[ix0:ix1], all_f[ix1:ix2], all_f[ix2::]


def check_mol2_line(file_name):
    with open(file_name, 'r') as f:
        while f.readline() != "@<TRIPOS>ATOM\n":
            pass
        return len(f.readline().split())


def xyz_names(file_name):
    positions = read_mol2(file_name)[1].split()[1::]
    names = {}
    xyz = {}
    ns = check_mol2_line(file_name)
    for i in range(len(positions)//ns):
        names.update({int(positions[ns*i]): positions[ns*i+1]})
        xyz.update({int(positions[ns*i]): np.array([float(positions[ns*i+2]),
                                                 float(positions[ns * i + 3]),
                                                 float(positions[ns * i + 4])])})
    return xyz, names



def xyz_names_bonds(file_name):
    _, positions, bondsf = read_mol2(file_name)
    positions = positions.split()[1::]
    bondsf = bondsf.split()[1::]
    names, xyz, bonds = {}, {}, {}
    ns = check_mol2_line(file_name)
    for i in range(len(positions) // ns):
        names.update({int(positions[ns * i]): positions[ns * i + 1]})
        xyz.update({int(positions[ns * i]): np.array([float(positions[ns * i + 2]),
                                                 float(positions[ns * i + 3]),
                                                 float(positions[ns * i + 4])])})
    for i in range(len(bondsf) // 4):
        b1, b2, attr = bondsf[4 * i + 1:4 * i + 4:]
    return xyz, names


if __name__ == "__main__":
    print(xyz_names_bonds('benzene.mol2'))
    print(check_mol2_line('benzene.mol2'))

