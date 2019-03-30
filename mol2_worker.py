import numpy as np
from copy import deepcopy
from quadro_with_rotate_class import Spherical_divider
import rmsd


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

    def set_xyz(self, x, y, z):
        self.z = z
        self.y = y
        self.x = x

    def defto_string(self):
        return '\t'.join(list(map(str, [self.name, self.x, self.y, self.z, self.name_i, self.i1, self.i2, self.i3])))


    def position(self):
        return np.array([self.x, self.y, self.z])


class Bond():
    def __init__(self, c1, c2, attr, length=1., sections=[]):
        self.connected = {c1, c2} # c1<c2
        self.attr = attr
        self.length = length
        self.sections = sections # [section_c1, section_c2]

    def set_length(self, length):
        self.length = length


    def set_section(self, sections):
        self.sections = sections


class Molecule():
    def __init__(self, mol2file='', n=5):
        self.atoms = {}
        self.n = n
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

    def set_n(self, n):
        self.n = n

    def set_notation(self):
        '''
        self.notation is array [section: from max num, section2: from min num]
        :return: None
        '''
        self.divider = Spherical_divider(self.n)
        for kbond, ibond in self.bonds.items():
            section = self.divider.find_section(self.atoms[min(ibond.connected)].position(), self.atoms[max(ibond.connected)].position())
            section_anti = self.divider.find_section(self.atoms[max(ibond.connected)].position(), self.atoms[min(ibond.connected)].position())
            ibond.set_section([section, section_anti])

    def check(self, dim_structure_reduced, eps=0.01, power_eps=3):
        forces = 1
        forces_next = 0
        bs = bonds_of_paired(self.bonds)
        # print(dim_structure_reduced)
        while abs(forces_next-forces) > eps ** power_eps:
            apply_force = {}
            for key, item in dim_structure_reduced.items():
                bonds_of_key_atom = [ar for _, ar in self.bonds.items() if (key in ar.connected
                       and list(ar.connected-{key})[0] in dim_structure_reduced.keys())]
                apply_force.update({key: np.array([0, 0, 0])})
                for b in bonds_of_key_atom:
                    another = list(b.connected-{key})[0]
                    apply_force[key] = apply_force[key] - dim_structure_reduced[key]+(dim_structure_reduced[another]
                                           +self.divider.scube[bs[tuple(sorted([key, another]))][1][key > another]])
                dim_structure_reduced[key] = dim_structure_reduced[key] + eps * apply_force[key]
                forces, forces_next = forces_next, sum([np.linalg.norm(el) for _, el in apply_force.items()])
        return dim_structure_reduced


    def from_notation(self):
        """
        From bonds cube notation get xyz-coords with relax
        :return dim_structure:
        """
        dim_structure = {1: np.array([0, 0, 0])}
        bonds = bonds_of_paired(self.bonds)
        bonds_copy = deepcopy(bonds)
        p = [[1, max(key), bonds_copy[key]] for key in bonds_copy.keys() if 1 in key]
        while len(p) != 0:
            cur_key, atom_to_bond, bond_characteristics = p.pop(0)
            if not (atom_to_bond in dim_structure.keys()):
                coord = self.divider.scube[bond_characteristics[1][cur_key < atom_to_bond]] * bond_characteristics[0]\
                        + dim_structure[cur_key]
                dim_structure.update({atom_to_bond: coord})
                [p.append([atom_to_bond, list(set(key)-{atom_to_bond})[0], bonds_copy[key]])
                 for key in bonds_copy.keys() if atom_to_bond in key]
            self.check(dim_structure, eps=0.1, power_eps=1)
        return dim_structure



    def to_mol2(self, mol2file):
        with open(mol2file, 'w') as f:
            f.write('@<TRIPOS>MOLECULE\n*****\n{} {} 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n'.format(str(len(self.atoms)), str(len(self.bonds))))
            for k, atom in self.atoms.items():
                f.write(('\t{}'*9+'\n').format(k, atom.name, atom.x, atom.y, atom.z, atom.name_i, atom.i1, atom.i2, atom.i3))
            f.write('@<TRIPOS>BOND\n')
            for k, bond in self.bonds.items():
                f.write('\t{}\t{}\t{}\t{}\n'.format(k, *bond.connected, bond.attr))


    def to_xyz(self, xyzfile):
        with open(xyzfile, 'w') as f:
            f.write(str(len(self.atoms))+'\n\n')
            for k, atom in self.atoms.items():
                f.write('{}\t{}\t{}\t{}\n'.format(atom.name, atom.x, atom.y, atom.z))

    def to_positions_array(self):
        pos = []
        for _, atom in self.atoms.items():
            pos.append(atom.position())
        return np.array(pos)


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

def bonds_of_paired(bonds):
    '''
    :param bonds: in format {1: [[[2,1]], (2,1)], ...}
    :return: bonds in format {(c_1, c_2): [length, sections, attr], ...}
    '''
    paired_bonds = {}
    for _, item in bonds.items():
        paired_bonds.update({tuple(sorted(list(item.connected))): [item.length, item.sections, item.attr]})
    return paired_bonds


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

def compare_structers(mol1, mol2):
    mol1 -= rmsd.centroid(mol1)
    mol2 -= rmsd.centroid(mol2)
    rotate = rmsd.kabsch(mol2, mol1)

    mol2 = np.dot(mol2, rotate)
    return rmsd.rmsd(mol1, mol2)

def xyz_to_array(file_name):
    with open(file_name, 'r') as f:
        n = int(f.readline())
        f.readline()
        lines = []
        for _ in range(n):
            lines.append([float(i) if len(i) > 3 and i[3:].isdigit() else i for i in f.readline().split()])
    return lines

if __name__ == "__main__":
    # bonds = (xyz_names_bonds('Caffein.mol2')[-1])
    # atoms = atoms_and_bonds('Caffein.mol2')

    a = Molecule('./many_opted_mol2s/3a-MnH2-ads-MeOAc_opted.mol2')

    a.set_notation()
    before = a.to_positions_array()
    struct = a.from_notation()

    b = deepcopy(a)
    for k, pos in struct.items():
        b.atoms[k].set_xyz(*pos)
    b.to_mol2('uncompressed.mol2')
    after = b.to_positions_array()

    print(compare_structers(before, after))

