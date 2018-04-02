import numpy as np

class Atom():
    def __init__(self, name, name_with_i, i1, i2, i3):
        self.name = name
        self.x = 0
        self.y = 0
        self.z = 0
        self.name_i = name_with_i
        self.i1 = i1
        self.i2 = i2
        self.i3 = i3
        self.orientation = (0, 0)

    def set_xyz(self, x, y, z):
        self.z = z
        self.y = y
        self.x = x

    def to_string(self):
        return '\t'.join(list(map(str, [self.name, self.x, self.y, self.z, self.name_i, self.i1, self.i2, self.i3])))

    def set_orientation(self, basis):
        self.orientation = basis

class Bond():
    def __init__(self, c1, c2):
        self.connected = {c1, c2}
        
def read_mol2(file_name):
    with open(file_name, 'r') as f1:
        all_f = f1.read()
    ix0 = all_f.index('@<TRIPOS>MOLECULE') # info = [ix0:ix1]
    ix1 = all_f.index('@<TRIPOS>ATOM')     # xyz_part = [ix1:ix2]
    ix2 = all_f.index('@<TRIPOS>BOND')     # bond = [ix2::]
    ix3 = all_f.find('@<TRIPOS>SUBSTRUCTURE')
    return all_f[ix0:ix1], all_f[ix1:ix2], all_f[ix2:ix3:]


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
    # print('b',bondsf)
    ns = check_mol2_line(file_name)
    for i in range(len(positions) // ns):
        names.update({int(positions[ns * i]): positions[ns * i + 1]})
        xyz.update({int(positions[ns * i]): np.array([float(positions[ns * i + 2]),
                                                 float(positions[ns * i + 3]),
                                                 float(positions[ns * i + 4])])})
    bonds = []
    for i in range(len(bondsf) // 4):
        b1, b2, attr = bondsf[4 * i + 1:4 * i + 4:]
        bonds.append([int(b1), int(b2), attr])
    return xyz, names, bonds


def positions_atoms_bonds(file_name):
    pass

if __name__ == "__main__":
    # print(xyz_names_bonds('benzene.mol2'))
    # print(check_mol2_line('benzene.mol2'))
    bonds = (xyz_names_bonds('Caffein.mol2')[-1])
    print(bonds)
    import pybel
