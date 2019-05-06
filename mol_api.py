import os
import copy
from tempfile import mkdtemp
from shutil import rmtree

from layouter import *
from mol2_chain_q import bonds_to_one_way_dict, to_two_ways_bond2, prepare_bonds
from mol2_worker import read_mol2, Atom, Bond, xyz_names_bonds, compare_structers
from mopac_worker import get_heat_of_xyz
from quadro_with_rotate_class import Spherical_divider


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
        for key in self.bonds.keys():
            self.bonds[key] = sorted(self.bonds[key])

    def diff(self, to_the_notation, b_diffs=False):
        '''
        TODO later it'd return the differents and check nums of atoms more correctly
        :param to_the_notation: compare with this one notation
        :param b_diffs: SOON
        :param s_diffs:
        :param l_diffs:
        :return: int of different bonds, int of different section, sum of different bonds lengths
        '''
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        sections_diffs = 0
        bonds_diffs = 0
        length_diff = 0
        for k_1, i_1 in self.notation.items():
            k_2 = k_1
            i_2 = to_the_notation.notation[k_2]
            if len(i_1[0]) == len(i_2[0]):
                for i, j in zip(i_1[0], i_2[0]):
                    if i[0] != j[0]:
                        bonds_diffs += 1
                        print(k_1, i,j)
                    elif i[1] != j[1]:
                        sections_diffs += 1
            else:
                bonds_diffs += abs(len(i_2)-len(i_1))
        for k, i in self.bonds.items():
            if to_the_notation.bonds.get(k):
                b2 = sorted(to_the_notation.bonds[k])
            else:
                continue
            i = sorted(i)
            if len(i) == len(b2):
                for b01, b02 in zip(i, b2):
                    if b01[0] == b02[0]:
                        length_diff += abs(b02[1]-b01[1])
            # else: calc smth else
        return bonds_diffs, sections_diffs, length_diff

    def s_diff(self, to_the_notation):
        '''
        :param to_the_notation: compare with this one notation
        :return:
        '''
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        s_d = []
        for k_1, i_1 in self.notation.items():
            k_2 = k_1
            i_2 = to_the_notation.notation[k_2]
            if len(i_1[0]) == len(i_2[0]):
                for inx, j in enumerate(zip(i_1[0], i_2[0])):
                    i, j = j
                    if i[0] == j[0] and i[1] != j[1]:
                        s_d.append([k_1, inx, j[1]])
        return s_d

    def l_diff(self, to_the_notation):
        '''
        :param to_the_notation: compare with this one notation
        :return:
        '''
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        l_d = []
        for k, i in self.bonds.items():
            if to_the_notation.bonds.get(k):
                b2 = sorted(to_the_notation.bonds[k])
            else:
                continue
            # i = sorted(i)
            if len(i) == len(b2):
                for inx, b01 in enumerate(zip(i, b2)):
                    b_1, b_2 = b01
                    if b_1[0] == b_2[0]:
                        l_d.append([k, inx, b_2])
            # else: calc smth else
        return l_d

    def get_heat(self, tmp=''):
        del_flag = False
        names = [i.name for _, i in self.atoms.items()]
        if tmp == '':
            tmp = mkdtemp()
            del_flag = True
        ds = [i for _, i in dimensional_structure(self).items()]
        file = os.path.join(tmp, 'to_calc_heat.xyz')
        with open(file, 'w') as f:
            n = len(names)
            f.write(str(n) + '\n\n')
            for ix in range(n):
                f.write('{}\t{}\t{}\t{}\n'.format(names[ix], *list(map(str, ds[ix]))))
        if del_flag: rmtree(tmp)
        return get_heat_of_xyz(file, tmpdir=tmp)


    def s_change(self, to_the_notation, follow_energy=False, sh=False, keep_change=False):
        if follow_energy:
            ens = []
            tmp = mkdtemp()
        s_d = self.s_diff(to_the_notation)
        if sh: np.random.shuffle(s_d)
        for k_1, inx, j in s_d:
            self.notation[k_1][0][inx][1] = j
            if follow_energy:
                ens.append(self.get_heat(tmp=tmp))
        if follow_energy:
            rmtree(tmp)
            return (ens, s_d) if keep_change else ens

    def s_change_step(self, to_the_notation):
        s_d = self.s_diff(to_the_notation)
        k_1, inx, j = s_d[0]
        self.notation[k_1][0][inx][1] = j

    def l_change_step(self, to_the_notation):
        l_d = self.l_diff(to_the_notation)
        k_1, inx, j = l_d[0] #TODO 0 will be a random index
        self.bonds[k_1][inx][1] = j[1]

    def l_change(self, to_the_notation, follow_energy=False, sh=False, keep_change=False):
        if follow_energy:
            ens = []
            tmp = mkdtemp()
        l_d = self.l_diff(to_the_notation)
        if sh: np.random.shuffle(l_d)
        for k_1, inx, j in l_d:
            self.bonds[k_1][inx][1] = j[1]
            if follow_energy:
                ens.append(self.get_heat(tmp=tmp))
        if follow_energy:
            rmtree(tmp)
            return (ens, l_d) if keep_change else ens


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

    def refresh_dimensional(self):
        ds = self.get_dimensional()
        for k, i in ds.items():
            self.atoms[k].x = i[0]
            self.atoms[k].y = i[1]
            self.atoms[k].z = i[2]

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
        :return: pair of atoms nums with zero bonds
        '''
        zero_pairs = []
        for k, i in self.bonds.items():
            if i.attr == '0':
                zero_pairs.append(i.connected)
        return zero_pairs

    def choose_bond(self, n=1):
        '''
        :return: pair of atoms nums to change bond
        '''
        bs = []
        for _ in range(n):
            r1 = np.random.choice([i for i in self.notation.bonds.keys()])
            r2 = np.random.choice([i[0] for i in self.notation.bonds[r1]])
            bs.append({r1, r2})
        return bs

    def get_index_of_bond(self, atom_with_bond, atom_to_bond):
        '''
        :param atom_with_bond: num of atom with bond
        :param atom_to_bond: num of atom bonded with atom_with_bond
        :return: index of atom_to_bond in notation[atom_with_bond]
        '''
        for inx, nbond in enumerate(self.notation.notation[atom_with_bond][0]):
            if nbond[0] == atom_to_bond:
                return inx
        return None

    def lazy_bond(self):
        '''
        :return: random atom and random index of bonded with it atom
        '''
        r1 = np.random.choice([i for i in self.notation.bonds.keys()])
        r2 = np.random.randint(len(self.notation.bonds[r1]))
        return r1, r2

    def get_child(self, zero=True, change_length=True, change_section=True, one_time_many=1):
        child = copy.deepcopy(self)
        if zero:
            bonds_to_change = self.zero_bonds()
        else:
            bonds_to_change = self.choose_bond(n=one_time_many)

        if change_section:
            for i, j in bonds_to_change:
                inx = self.get_index_of_bond(i, j)
                n_section = np.random.randint(0, len(self.notation.divider.scube))
                na_section = self.notation.divider.anti_scube[n_section]
                child.notation.notation[i][0][inx][1] = n_section
                inx1 = self.get_index_of_bond(j, i)
                child.notation.notation[j][0][inx1][1] = na_section

        if change_length: # it has collisions
            if zero:
                pass
            else:
                r1, r2 = self.lazy_bond()
                child.notation.bonds[r1][r2][1] += -0.1 if np.random.randint(2) else 0.1
        return child # change both length but how????

def searcher(substrat, product,zero_bonds_only=False,
             length_change=True, length_path=10):
    '''
    :param substrat: type: Molecule,
    :param product: type: Molecule,
    :param t_times: allowed rmsd increase
    :return: path of rmsd
    '''
    # c1 = substrat.compare_with(product.to_positions_array())
    dim_structure = dimensional_structure(substrat.notation, relax=True)
    c1 = product.compare_with(np.array([i for _, i in dim_structure.items()]))
    paths = [c1]
    st = copy.deepcopy(substrat)
    s_comp = substrat.notation.diff(product.notation)[1]
    while len(paths) < length_path:
    # while paths[-1] > 0.25:

        ch = st.get_child(zero=zero_bonds_only, change_length=length_change)
        ch_dim = dimensional_structure(ch.notation, relax=True)
        compared = pr.compare_with(np.array([i for _, i in ch_dim.items()]))
        _, s_comp_c, _ = st.notation.diff(pr.notation)
        if s_comp_c < s_comp:
        # if 1-np.random.rand() < c1/compared:
            paths.append(compared)
            print('structure is ', ch_dim)
            del st
            st = copy.deepcopy(ch)
    return paths

def random_search(reactant, product):
    '''
    :param reactant:
    :param product:
    :return: energies_path and length of sections changes
    '''
    s = reactant.notation.s_change(product.notation, follow_energy=True, sh=True)
    l = reactant.notation.l_change(product.notation, follow_energy=True, sh=True)
    return s+l, len(s)


if __name__ == '__main__':
    # params = {'n': 1,
    #           'reaction': 'mopac_example', #'3a->4'
    #           }
    n = 10
    reaction = 'mopac_example' # '3a->4' #
    # reaction = '3a->4' #
    if reaction == '3a->4':
        ln = Molecule('./ordered_mol2/3a.mol2', n=n)
        pr = Molecule('./ordered_mol2/4_opted.mol2', n=n)
    else:
        ln = Molecule('./ordered_mol2/js_exapmle_init.mol2', n=n)
        pr = Molecule('./ordered_mol2/js_exapmle_finish.mol2', n=n)
    # path = searcher(ln, pr)
    # print(path)
    # print(ln.notation.diff(pr.notation))
    # print('reactant with orignal reactant')
    # print(ln.compare_with(np.array([i for _, i in dimensional_structure(ln.notation).items()])))
    # print('before original pr with orignal reactant')
    # print(compare_structers(ln.to_positions_array(), pr.to_positions_array()))
    # print('before original pr with dim_struct of reactant')
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(ln.notation).items()])))
    # print('before dim_structure of pr with dim_struct of reactant')
    # print(compare_structers(np.array([i for _, i in dimensional_structure(ln.notation).items()]), np.array([i for _, i in dimensional_structure(pr.notation).items()])))
    # print('product with product')
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(pr.notation).items()])))
    # name_heat = reaction + 'start.xyz'
    # pr.to_xyz(name_heat)
    # print('before', get_heat_of_xyz(name_heat))
    # pr.refresh_dimensional()
    # name_heat = reaction+'start.xyz'
    # pr.to_xyz(name_heat)
    # print('after', get_heat_of_xyz(name_heat))
    print('length+sections')
    print(ln.notation.s_change(pr.notation, follow_energy=True, sh=True))
    # print(ln.notation.l_change(pr.notation, follow_energy=True))
    # print(ln.notation.s_change(pr.notation, follow_energy=True))

    # n = 10
    # reaction = 'mopac_example'  # '3a->4' #
    # # reaction = '3a->4' #
    # if reaction == '3a->4':
    #     ln = Molecule('./ordered_mol2/3a.mol2', n=n)
    #     pr = Molecule('./ordered_mol2/4_opted.mol2', n=n)
    # else:
    #     ln = Molecule('./ordered_mol2/js_exapmle_init.mol2', n=n)
    #     pr = Molecule('./ordered_mol2/js_exapmle_finish.mol2', n=n)
    # print('sections+length')
    # print(ln.notation.s_change(pr.notation, follow_energy=True))
    # print(ln.notation.l_change(pr.notation, follow_energy=True))

    # print(ln.notation.s_diff(pr.notation))
    # print(ln.notation.l_diff(pr.notation))

    # for k, i in pr.notation.bonds.items():
    #     print('atom ', k)
    #     print('product', i)
    #     print('reactant', ln.notation.bonds[k])
    # print(pr.notation.bonds)
    # print(ln.notation.l_diff(pr.notation))
    # print(dimensional_structure(ln.notation))

    # print(compare_structers(np.array([i for _, i in dimensional_structure(ln.notation).items()]), np.array([i for _, i in dimensional_structure(pr.notation).items()])))
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(ln.notation).items()])))
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(ln.notation).items()])))
    # print(ln.compare_with(pr.to_positions_array()))


    # print(ln.compare_with(pr.to_positions_array()))
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(ln.notation).items()])))

    # plt.plot(list(range(len(paths))), paths)
    # plt.title('reaction {},n={}, all changes'.format(reaction, str(n)))
    # plt.show()
    # write_mol2_file('My_' + name + '_' + 'q0.mol2', atoms_info, dim_structure, bonds=paired)