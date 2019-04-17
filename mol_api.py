import numpy as np
from quadro_with_rotate_class import Spherical_divider
from mol2_worker import bonds_of_paired, read_mol2, Atom, Bond, xyz_names_bonds, compare_structers
from mol2_chain_q import bonds_to_one_way_dict, to_two_ways_bond2, prepare_bonds
from layouter import *
import copy

class Notation():
    def __init__(self, n, info_from_file):
        self.divider = Spherical_divider(n=n)
        bonds, self.atoms = info_from_file
        positions_copy = copy.deepcopy(self.atoms)
        self.notation = {}
        self.bonds = bonds_to_one_way_dict(prepare_bonds(bonds))
        bonds2 = to_two_ways_bond2(self.bonds, with_attr=True)
        for key, item in self.atoms.items():
            cur_p = positions_copy.pop(key).position()
            connected = [i[0] for i in bonds2[key]]
            self.notation.update({key: [list([i, self.divider.find_section(cur_p, self.atoms[i].position())]
                                        for i in connected)]})
        for key, item in self.bonds.items():
            for i in range(len(item)):
                self.bonds[key][i].insert(1, round(np.linalg.norm(self.atoms[key].position() - self.atoms[item[i][0]].position()), 1))


class Molecule():
    def __init__(self, mol2file='', n=5):
        self.atoms = {}
        self.n = n
        self.mol2file = mol2file
        info, positions, bonds = read_mol2(self.mol2file)
        for line in positions.split('\n')[1:-1:]:
            l = line.split()
            self.atoms.update({int(l[0]): Atom(l[1], *l[5::])})
            self.atoms[int(l[0])].set_xyz(*list(map(float, l[2:5:])))
        self.notation = Notation(n=self.n, info_from_file=xyz_names_bonds(self.mol2file))
        self.bonds = {}
        for line in bonds.split('\n')[1::]:
            l = line.split()
            l = list(map(int, l[:3:]))+[l[3]]
            self.bonds.update({l[0]: Bond(*l[1::],
                                          length=np.linalg.norm(self.atoms[l[1]].position()
                                                                - self.atoms[l[2]].position()))})
    def set_n(self, n):
        self.n = n
        self.notation = Notation(n=n, info_from_file=xyz_names_bonds(self.mol2file))

    def get_dimensional(self, relax=True):
        return dimensional_structure(self.notation, relax=relax)


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

    def compare_with(self, xyz_positions):
        origin = copy.deepcopy(self.to_positions_array())
        xyz_positions_c = copy.deepcopy(xyz_positions)
        return compare_structers(origin, xyz_positions_c)

    def zero_bonds(self):
        '''
        :return: pair of atoms with zero bonds
        '''
        zero_pairs = []
        for k, i in self.bonds.items():
            if i.attr == '0':
                zero_pairs.append(i.connected)
        return zero_pairs

    def get_child(self, zero=True):
        zero_bonds = self.zero_bonds()
        child = copy.deepcopy(self)
        for i, j in zero_bonds:
            for inx, nbond in enumerate(self.notation.notation[i][0]):
                if nbond[0] == j:
                    break
            n_section = np.random.randint(0, len(self.notation.divider.scube))
            na_section = self.notation.divider.anti_scube[n_section]
            child.notation.notation[i][0][inx][1] = n_section
            for inx, nbond in enumerate(self.notation.notation[j][0]):
                if nbond[0] == i:
                    break
            child.notation.notation[j][0][inx][1] = na_section

        ### and change length of bond ### fix it

        # r1 = np.random.randint(1, len(self.notation.bonds.keys())+1)
        # r2 = np.random.randint(len(self.notation.bonds[r1]))
        # child.notation.bonds[r1][r2][1] += -0.1 if np.random.randint(1) else 0.1

        return child # change both length but how????


if __name__ == '__main__':
    ln = Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/js_exapmle_init.mol2', n=3)
    pr = Molecule('/home/anastasiia/PycharmProjects/chem_structure/ordered_mol2/js_exapmle_finish.mol2', n=3)
    print('init to finish rmds ', ln.compare_with(pr.to_positions_array()))
    # print(ln.notation.notation)
    # ln = Notation(n=n_param, info_from_file=xyz_names_bonds(file_name + '.mol2'))
    # paired = bonds_of_paired(ln.bonds)
    dim_structure = dimensional_structure(ln.notation, relax=True)
    # print(ln.to_positions_array())
    # print(ln.compare_with(np.array([i for _, i in dim_structure.items()])))

    ch = ln.get_child()
    ch_dim = dimensional_structure(ch.notation, relax=True)
    print('child to finish rmds ', pr.compare_with(np.array([i for _, i in ch_dim.items()])))

    # for _ in range(10):

    print(ln.compare_with(np.array([i for _, i in ch_dim.items()])))

    # write_mol2_file('My_' + name + '_' + 'q0.mol2', atoms_info, dim_structure, bonds=paired)