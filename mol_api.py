import os
from tempfile import mkdtemp
from shutil import rmtree
import matplotlib.pyplot as plt

from layouter import *
from mol2_worker import read_mol2, Atom, xyz_names_bonds, compare_structers, Bond, bonds_to_dict
from mopac_worker import get_heat_of_xyz, get_energy_of_xyz
from quadro_with_rotate_class import Spherical_divider


class Notation:
    def __init__(self, n, info_from_file):
        self.divider = Spherical_divider(n=n)
        bonds, self.atoms = info_from_file
        self.notation = {}
        self.bonds = bonds_to_dict(bonds)
        self.set_notation(copy.deepcopy(self.atoms))

    def set_notation(self, positions_copy):
        for key, item in self.atoms.items():
            cur_p = positions_copy.pop(key).position()
            self.notation.update({key: {}})
            for k, _ in self.bonds[key].items():
                self.notation[key].update({k: self.divider.find_section(cur_p, self.atoms[k].position())})
                if self.notation.get(k) is None:
                    self.notation.update({k: {}})
                self.notation[k].update({key: self.divider.find_section(self.atoms[k].position(), cur_p)})
                l = round(np.linalg.norm(self.atoms[key].position() - self.atoms[k].position()), 1)
                self.bonds[key][k].set_length(l)
                self.bonds[key][k].set_section(self.notation[key][k])

    def difference_of_bonds(self, to_the_notation):
        bonds_diffs = 0
        for k_1, i_1 in self.notation.items():
            if to_the_notation.notation.get(k_1) is None:
                bonds_diffs += len(i_1.keys())
            else:
                i_2 = set(to_the_notation.notation.get(k_1).keys())
                bonds_diffs += len(set(i_1.keys()) ^ i_2)
        d = set(to_the_notation.notation.keys()) - set(self.notation.keys())
        bonds_diffs += sum([len(to_the_notation[dd]) for dd in d])
        return bonds_diffs

    def difference_of_sections(self, to_the_notation):
        sections_diffs = 0
        for k_1, i_1 in self.notation.items():
            if to_the_notation.notation.get(k_1) is None:
                continue
            k_2 = k_1
            ki1 = list(i_1.keys())
            i_2 = to_the_notation.notation.get(k_2)
            s1 = np.array([i for _, i in i_1.items()])
            s2 = np.array([i_2.get(inx) for inx in ki1])
            sections_diffs += sum(np.where(s1 == s2, 0, 1))
        return sections_diffs

    def difference_of_lengths(self, to_the_notation):
        length_diff = 0
        for k, i in self.bonds.items():
            if not to_the_notation.bonds.get(k):
                continue
            for k2, i2 in self.bonds[k].items(): # k - k2 pairs
                if to_the_notation.bonds[k].get(k2):
                    length_diff += abs(i2.length - to_the_notation.bonds[k][k2].length)
        return length_diff

    def diff(self, to_the_notation):
        """
        :param to_the_notation: compare with this one notation
        :return: int of different bonds, int of different section, sum of different bonds lengths
        """
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        return self.difference_of_bonds(to_the_notation), self.difference_of_sections(to_the_notation),\
               self.difference_of_lengths(to_the_notation)

    def s_diff(self, to_the_notation):
        """
        :param to_the_notation: compare with this one notation
        :return:
        """
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        s_d = []
        for k_1, i_1 in self.notation.items():
            i_2 = to_the_notation.notation.get(k_1)
            if i_2 is None:
                return []
            ki1 = list(i_1.keys())
            s1 = np.array([i for _, i in i_1.items()])
            s2 = np.array([i_2.get(inx) for inx in ki1])
            for a, sec1, sec2 in zip(i_1.items(), s1, s2):
                if (not sec2 is None) and sec1 != sec2:
                    s_d.append([k_1, a[0], sec2])
        return s_d

    def l_diff(self, to_the_notation):
        """
        :param to_the_notation: compare with this one notation
        :return:
        """
        if self.divider.n != to_the_notation.divider.n:
            print("Different dividers!")
            return
        l_d = []
        for k, i in self.bonds.items():
            if to_the_notation.bonds.get(k):
                b2 = to_the_notation.bonds[k]
            else:
                continue
            for k2, i2 in i.items():
                if b2.get(k2):
                    to_l = b2[k2].length
                    if i2.length != to_l: l_d.append([k, k2, to_l])
        return l_d

    def get_heat(self, tmp=''):
        del_flag = False
        names = [i.name for _, i in self.atoms.items()]
        if tmp == '':
            tmp = mkdtemp()
            del_flag = True
        file = os.path.join(tmp, 'to_calc_heat.xyz')
        with open(file, 'w') as f:
            n = len(names)
            f.write(str(n) + '\n\n')
            for ix in range(n):
                f.write('{}\t{}\t{}\t{}\n'.format(names[ix], self.atoms[ix+1].x, self.atoms[ix+1].y, self.atoms[ix+1].z))
        heat = get_heat_of_xyz(file, tmpdir=tmp)
        if del_flag: rmtree(tmp)
        return heat


    def s_change(self, to_the_notation, follow_energy=False, sh=False, keep_change=False):
        if follow_energy:
            ens = []
            tmp = mkdtemp()
        s_d = self.s_diff(to_the_notation)
        if sh: np.random.shuffle(s_d)
        for k_1, inx, j in s_d:
            self.notation[k_1][inx] = j
            if follow_energy:
                ens.append(self.get_heat(tmp=tmp))
        if follow_energy:
            rmtree(tmp)
            return (ens, s_d) if keep_change else ens

    def s_change_step(self, to_the_notation):
        s_d = self.s_diff(to_the_notation)
        if s_d != []:
            indx = np.random.randint(len(s_d))
            k_1, inx, j = s_d[indx]
            self.notation[k_1][inx] = j
            self.bonds[k_1][inx].set_section(j)
            self.bonds[inx][k_1].set_section(self.divider.anti_scube[j])
            return 0
        return -1

    def l_change_step(self, to_the_notation):
        l_d = self.l_diff(to_the_notation)
        if l_d != []:
            indx = np.random.randint(len(l_d))
            k_1, inx, j = l_d[indx]
            self.bonds[k_1][inx].set_length(j)
            return 0
        return -1

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


class Molecule:
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
        self.bonds_from_mol2 = {}
        for line in bonds.split('\n')[1::]:
            l = line.split()
            l = list(map(int, l[:3:]))+[l[3]]
            self.bonds_from_mol2.update({l[0]: Bond(*l[1::],
                                          length=np.linalg.norm(self.atoms[l[1]].position()
                                                                - self.atoms[l[2]].position()))})

    def refresh_dimensional(self, relax=True):
        ds = self.get_dimensional(relax=relax)
        for k, i in ds.items():
            self.atoms[k].x = i[0]
            self.atoms[k].y = i[1]
            self.atoms[k].z = i[2]

    def set_dimensional(self, positions):
        for k, i in enumerate(positions):
            self.atoms[k+1].x = i[0]
            self.atoms[k+1].y = i[1]
            self.atoms[k+1].z = i[2]

    def set_n(self, n):
        self.n = n
        self.notation = Notation(n=n, info_from_file=xyz_names_bonds(self.mol2file))

    def add_bond(self, c1, c2, length, section, attr='1'): #TODO fix this one function
        '''
        Warning: use only if unatural add (otherwise calculate both sections)
        :param c1, c2: c1 < c2
        :param attr: string_describe
        '''
        asection = self.notation.divider.anti_scube[section]
        if self.notation.notation.get(c1) == None:
            self.notation.notation.update({c1: {}})
        self.notation.notation[c1].update({c2: section})
        if self.notation.notation.get(c2) == None:
            self.notation.notation.update({c2: {}})
        self.notation.notation[c2].update({c1: asection})
        if self.notation.bonds.get(c1) == None:
            self.notation.bonds.update({c1: {}})
        if self.notation.bonds.get(c2) == None:
            self.notation.bonds.update({c2: {}})
        self.notation.bonds[c1].update({c2: Bond(c1, c2, attr=attr, length=length, section=section)})
        self.notation.bonds[c2].update({c1: Bond(c2, c1, attr=attr, length=length, section=asection)})
        self.notation.notation[c1].update({c2: section})
        self.notation.notation[c2].update({c1: asection})

    def check_c1_section_is_free(self, c1, section):
        sections = [i for _, i in c1.items()]
        return not (section in sections)

    def get_dimensional(self, relax=True):
        return dimensional_structure(self.notation, relax=relax)

    def to_mol2(self, mol2file):
        with open(mol2file, 'w') as f:
            f.write('@<TRIPOS>MOLECULE\n*****\n{} {} 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n'.format(str(len(self.atoms)), str(len(self.bonds_from_mol2))))
            for k, atom in self.atoms.items():
                f.write(('\t{}'*9+'\n').format(k, atom.name, atom.x, atom.y, atom.z, atom.name_i, atom.i1, atom.i2, atom.i3))
            f.write('@<TRIPOS>BOND\n')
            for k, bond in self.bonds_from_mol2.items():
                f.write('\t{}\t{}\t{}\t{}\n'.format(k, bond.c1, bond.c2, bond.attr))

    def to_xyz(self, xyzfile, title='', mode='w'):
        with open(xyzfile, mode) as f:
            f.write(str(len(self.atoms))+'\n{}\n'.format(title))
            for k, atom in self.atoms.items():
                f.write('{}\t{}\t{}\t{}\n'.format(atom.name, atom.x, atom.y, atom.z))

    def to_positions_array(self):
        pos = []
        for _, atom in self.atoms.items():
            pos.append(atom.position())
        return np.array(pos)

    def get_energy(self, tmp=''):
        del_flag = False
        names = [i.name for _, i in self.atoms.items()]
        if tmp == '':
            tmp = mkdtemp()
            del_flag = True
        file = os.path.join(tmp, 'to_calc_heat.xyz')
        with open(file, 'w') as f:
            n = len(names)
            f.write(str(n) + '\n\n')
            for ix in range(n):
                f.write('{}\t{}\t{}\t{}\n'.format(names[ix], self.atoms[ix+1].x, self.atoms[ix+1].y, self.atoms[ix+1].z))
        heat = get_energy_of_xyz(file, tmpdir=tmp)
        if del_flag: rmtree(tmp)
        return heat

    def interact_pair(self, distance=1.5, no_less_than=0.5):
        # TODO
        pass

    def compare_with(self, xyz_positions):
        origin = copy.deepcopy(self.to_positions_array())
        xyz_positions_c = copy.deepcopy(xyz_positions)
        return compare_structers(origin, xyz_positions_c)

    def zero_bonds(self):
        '''
        :return: pair of atoms nums with zero bonds
        '''
        zero_pairs = []
        for _, i in self.bonds.items():
            for _, iinner in i.items():
                if iinner.attr == '0':
                    zero_pairs.append(i.connected)
        return zero_pairs

    def choose_bond(self, n=1):
        '''
        :return: pair of atoms nums to change bond
        '''
        bs = []
        for _ in range(n):
            r1 = np.random.choice(self.notation.notation.keys())
            r2 = np.random.choice(self.notation.notation[r1].keys())
            bs.append({r1, r2})
        return bs

    def atom_and_bonded_with_it(self):
        '''
        :return: random atom and random index of bonded with it atom
        '''
        r1 = np.random.choice(list(self.notation.notation.keys()))
        r2 = np.random.choice(list(self.notation.notation[r1].keys()))
        return r1, r2

    def get_child(self, zero=True, change_length=True,
                  change_section=True, add_or_remove=True, one_time_many=1):
        child = copy.deepcopy(self)
        if zero: #TODO fix zero if it needs...
            bonds_to_change = self.zero_bonds()
        else:
            bonds_to_change = self.choose_bond(n=one_time_many)

        if change_section:
            for i, j in bonds_to_change:
                n_section = np.random.randint(0, len(self.notation.divider.scube))
                na_section = self.notation.divider.anti_scube[n_section]
                child.notation.notation[i][j] = n_section
                child.notation.notation[j][i] = na_section

        if change_length: # it has collisions
            if zero:
                pass
            else:
                r1, r2 = self.atom_and_bonded_with_it()
                child.notation.bonds[r1][r2].length = child.notation.bonds[r1][r2].length (-0.1 if np.random.randint(2) else 0.1)
        if add_or_remove:
            a1, a2 = np.random.choice(range(len(self.atoms)), 2, replace=False) + np.array([1,1])
            self.add_bond(a1, a2, round(np.random.normal(1.1, 0.3), 2), np.random.randint(len(self.notation.divider.scube)))
        return child

    def mutation(self):
        mutant = copy.deepcopy(self)
        a1, a2 = self.atom_and_bonded_with_it()
        choi = np.random.random()
        if choi < 0.4:
            n_section = np.random.randint(0, len(self.notation.divider.scube))
            while not self.check_c1_section_is_free(self.notation.notation[a1], n_section):
                n_section = np.random.randint(0, len(self.notation.divider.scube))
            mutant.notation.notation[a1][a2] = n_section
            mutant.notation.notation[a2][a1] = self.notation.divider.anti_scube[n_section]
        elif choi < 0.8:
            dl = np.random.normal(0, 0.5, 1)[0]
            mutant.notation.bonds[a1][a2].set_length(mutant.notation.bonds[a1][a2].length + round(dl, 1))
        else:
            b1, b2 = sorted(np.random.choice(range(len(self.atoms)), 2) + np.array([1, 1]))
            mutant.add_bond(b1, b2, round(np.random.normal(1.1, 0.3), 1),
                          np.random.randint(len(mutant.notation.divider.scube)))
        return mutant


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

def real_random_path(reactant, product, n=10000, write=False, file_log='random_search_report'):
    '''
    :param reactant, product: Molecule
    :param n: number of loops
    :return: show plot of rmsds
    '''
    def apply_change():
        mut.refresh_dimensional()
        # path.append(reactant.notation.get_heat())
        msds.append(compare_structers(mut.to_positions_array(), product.to_positions_array()))
        if write:
            mut.to_xyz(file_log, mode='a')
            with open(file_log, 'a') as f_w:
                f_w.write('################\n')
    with open(file_log, 'w') as f_w:
        f_w.write('\n\n\n\n')
    mut = copy.deepcopy(reactant)
    d = reactant.notation.diff(product.notation)
    msd = compare_structers(mut.to_positions_array(), product.to_positions_array())
    mut_pr = reactant.mutation()
    mut_pr.refresh_dimensional()
    d_pr = mut_pr.notation.diff(product.notation)
    msd_pr = compare_structers(mut_pr.to_positions_array(), product.to_positions_array())
    msds = [msd]
    apply_change()
    print('initial msd', msd, 'd', d)
    print('rmsd accept probability')
    for _ in range(n):
        # if d_pr[1] < d[1] or d_pr[2] < d[2] or np.random.random() < np.exp(-(d_pr[2]/d[2])):
        if d_pr[1] < d[1] or d_pr[2] < d[2] or msd_pr < msds[0]*1.5:
            mut, d, msd = copy.deepcopy(mut_pr), d_pr, msd_pr
            msds.append(msd)
            apply_change()
        mut_pr = mut.mutation()
        mut_pr.refresh_dimensional()
        d_pr = mut_pr.notation.diff(product.notation)
        msd_pr = compare_structers(mut_pr.to_positions_array(), product.to_positions_array())
    print('final msd', msd, 'd', d)
    if write:
        f = open(file_log, 'r+')
        f.seek(0, 0)
        f.write(str(len(msds)-1)+'\n')
        f.close()
    plt.plot(list(range(len(msds))), msds)
    plt.show()


def random_to_the_aim_search(reactant, product, write=False, file_log='reaction_report'): #TODO implement here add of bonds
    """
    :return: energies_path and length of sections changes
    """
    def apply_change():
        reactant.refresh_dimensional()
        # path.append(reactant.notation.get_heat())
        msds.append(compare_structers(reactant.to_positions_array(), product.to_positions_array()))
        if write:
            reactant.to_xyz(file_log, mode='a')
            with open(file_log, 'a') as f_w:
                f_w.write('################\n')
    with open(file_log, 'w') as f_w:
        f_w.write('153\n')
    d = reactant.notation.diff(product.notation)
    # path = [reactant.notation.get_heat()]
    msds = []
    apply_change()
    while d[2] > 0.1 and d[1] != 0:
        if np.random.random() < 0.5:
            reactant.notation.s_change_step(product.notation)
        else:
            reactant.notation.l_change_step(product.notation)
        apply_change()
        d = reactant.notation.diff(product.notation)
    while reactant.notation.l_change_step(product.notation) != -1:
        apply_change()
    while reactant.notation.s_change_step(product.notation) != -1:
        apply_change()
    if False and write:
        f = open(file_log, 'r+')
        f.seek(0, 0)
        f.write(str(len(msds))+'\n')
        f.close()

    # return path, msds
    return msds

def genetic_to_the_aim(reactant, product, write=False,
                             file_log='reaction_report'):
    """
    :return: energies_path and length of sections changes
    """
    def apply_change():
        mutant.refresh_dimensional()
        delta = mutant.get_energy() - path[-1]
        if np.random.random() < np.exp(-delta/20.):
            path.append(path[-1]+delta)
            if write:
                mutant.to_xyz(file_log, title=str(path[-1]), mode='a')
                with open(file_log, 'a') as f_w:
                    f_w.write('################\n')
            return True
        return False

    d = reactant.notation.diff(product.notation)
    with open(file_log, 'w') as f_w:
        f_w.write(str(len(d))+'\n')
    path = [reactant.get_energy()]
    # msds = []
    # apply_change()
    while d[2] > 0.1 and d[1] != 0:
        mutant = copy.deepcopy(reactant)
        if np.random.random() < 0.5:
            mutant.notation.s_change_step(product.notation)
        else:
            mutant.notation.l_change_step(product.notation)
        if apply_change(): reactant = copy.deepcopy(mutant)
        d = reactant.notation.diff(product.notation)
    while reactant.notation.l_change_step(product.notation) != -1:
        mutant = copy.deepcopy(reactant)
        if apply_change(): reactant = copy.deepcopy(mutant)
    while reactant.notation.s_change_step(product.notation) != -1:
        mutant = copy.deepcopy(reactant)
        if apply_change(): reactant = copy.deepcopy(mutant)
    return path #, msds
    # return msds

    #########two halves of random path#########
    # s = reactant.notation.s_change(product.notation, follow_energy=True, sh=True)
    # l = reactant.notation.l_change(product.notation, follow_energy=True, sh=True)
    # return s+l, len(s)
    ###########################################



if __name__ == '__main__':
    n = 15
    # reaction = 'mopac_example' # '3a->4' #
    reaction = '3a->4' #
    # reaction = 'vanadii'
    if reaction == '3a->4':
        ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n)
        pr = Molecule('./prepared_mols2/4_opted.mol2', n=n)
        pr.refresh_dimensional()
    elif reaction == 'vanadii':
        ln = Molecule('./vanadii/3a_singlet_opted.mol2', n=n)
        pr = Molecule('./vanadii/ts_3a_4a_opted.mol2', n=n)
        pr.refresh_dimensional()
    else:
        ln = Molecule('./ordered_mol2/js_exapmle_init.mol2', n=n)
        pr = Molecule('./ordered_mol2/js_exapmle_finish.mol2', n=n)
        pr.refresh_dimensional()
    # ll = ln.to_positions_array()
    # pp = pr.to_positions_array()
    # import rmsd
    # ll -= rmsd.centroid(ll)
    # pp -= rmsd.centroid(pp)
    # rotate = rmsd.kabsch(pp, ll)
    # pp = np.dot(pp, rotate)
    # ln.set_dimensional(ll)
    # pr.set_dimensional(pp)
    # ln.to_mol2('./vanadii/3a_singlet_rot.mol2')
    # pr.to_mol2('./vanadii/ts_3a_4_rot.mol2')
    # compare_structers()
    # print(ln.notation.get())
    maxs = []
    for i in range(50, 81):
        ln = Molecule('./prepared_mols2/3a_opted.mol2', n=n)
        pr = Molecule('./prepared_mols2/4_opted.mol2', n=n)
        pr.refresh_dimensional()
        ms = genetic_to_the_aim(ln, pr, write=True, file_log=reaction+'_gen_report_20_'+str(i))
        # print(ms, )
        maxs.append(max(ms))
    print(maxs)
    # print(ln.notation.get_energy())
    # print(pr.notation.get_energy())
    # real_random_path(ln, pr, n=10000, write=True, file_log=reaction+'_random_report')
    # print(max(p))
    # print(p)
    # print(max(ms))
    # print(ms)
    # print(pr.compare_with(np.array([i for _, i in dimensional_structure(pr.notation).items()])))
    # name_heat = reaction + 'start.xyz'
    # pr.to_xyz(name_heat)
    # print('before', get_heat_of_xyz(name_heat))
    # pr.refresh_dimensional()
    # name_heat = reaction+'start.xyz'
    # pr.to_xyz(name_heat)
    # print(ln.notation.l_change(pr.notation, follow_energy=True))
    # print(ln.notation.s_change(pr.notation, follow_energy=True))


    # write_mol2_file('My_' + name + '_' + 'q0.mol2', atoms_info, dim_structure, bonds=paired)