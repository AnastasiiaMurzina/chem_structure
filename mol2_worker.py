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

    def defto_string(self):
        return '\t'.join(list(map(str, [self.name, self.x, self.y, self.z, self.name_i, self.i1, self.i2, self.i3])))

    def set_orientation(self, basis):
        self.orientation = basis

    def position(self):
        return np.array([self.x, self.y, self.z])


class Bond():
    def __init__(self, c1, c2, attr, length=1., sections=0):
        self.connected = {c1, c2} # c1<c2
        self.attr = attr
        self.length = length
        self.sections = sections # [section_c1, section_c2]

    def set_length(self, length):
        self.length = length


    def set_section(self, sections):
        self.sections = sections


class Molecule():
    def __init__(self, mol2file=''):
        self.atoms = {}
        info, positions, bonds = read_mol2(mol2file)
        for line in positions.split('\n')[1:-1:]:
            l = line.split()
            self.atoms.update({int(l[0]): Atom(l[1], *l[5::])})
            self.atoms[int(l[0])].set_xyz(*list(map(float, l[2:5:])))
        self.bonds = {}
        for line in bonds.split('\n')[1::]:
            l = line.split()
            l = [int(l[0]), int(l[1]), int(l[2]), l[3]]
            self.bonds.update({l[0]: Bond(l[1], l[2], l[3])})
            self.bonds[l[0]].set_length(np.linalg.norm(self.atoms[l[1]].position()-self.atoms[l[2]].position()))

    def Notation(self):
        pass




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
    return xyz # , names


def xyz_names_bonds(file_name):
    _, positions, bondsf = read_mol2(file_name)
    positions = positions.split()[1::]
    bondsf = bondsf.split()[1::]
    atoms, xyz, bonds = {}, {}, {}
    ns = check_mol2_line(file_name)
    for i in range(len(positions) // ns):
        atom = Atom(positions[ns * i + 1], positions[ns*i + 5],
                    positions[ns*i + 6], positions[ns*i + 7],
                    float(positions[ns*i + 8]))
        atom.set_xyz(float(positions[ns * i + 2]), float(positions[ns * i + 3]), float(positions[ns * i + 4]))
        atoms.update({i+1: atom})
    bonds = []
    for i in range(len(bondsf) // 4):
        b1, b2, attr = bondsf[4 * i + 1:4 * i + 4:]
        bonds.append([int(b1), int(b2), attr])
    return bonds, atoms


def atoms_and_bonds(file_name, bonds_choice=False):
    atoms = {}
    bonds = {}
    _, atomsxyz, bonds_attr = read_mol2(file_name)
    atomsxyz = atomsxyz.split()[1::]
    ago = check_mol2_line(file_name)
    n = int(atomsxyz[-ago])
    for i in range(n):
        at = atomsxyz[ago*i:ago*(i+1):]
        current_atom = Atom(at[1], at[5], at[6], at[7], float(at[8]))
        current_atom.set_xyz(float(at[2]), float(at[3]), float(at[4]))
        atoms.update({i+1: current_atom})
    if bonds_choice:
        bonds_attr = bonds_attr.split()[1::]
        for i in range(len(bonds_attr)//4):
            bond = bonds_attr[4*i: 4*i + 4:]
            bonds.update({i+1: Bond(int(bond[1]), int(bond[2]), bond[3])})
        return atoms, bonds
    return atoms

if __name__ == "__main__":
    # bonds = (xyz_names_bonds('Caffein.mol2')[-1])
    # atoms = atoms_and_bonds('Caffein.mol2')
    a = Molecule('./many_opted_mol2s/3a-MnH2-ads-MeOAc_opted.mol2')
    print(a.bonds[1].connected)