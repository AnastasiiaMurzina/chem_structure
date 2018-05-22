from penta_with_rotate import get_penta_points, find_section,\
    rotate_by_basis, find_basis, rotate_non_perpendicular, rotate_ten_vars
from mol2_worker import xyz_names, xyz_names_bonds, Atom, atoms_and_bonds, Bond
import numpy as np
import copy

eps_length = 0.001
pp = get_penta_points()

###########################Builder structure if possile else return fail########################
def prepare_bonds(bonds):
    for i, j in enumerate(bonds):
        if j[0] > j[1]:
            bonds[i][0], bonds[i][1] = bonds[i][1], bonds[i][0]
    return bonds


def bonds_of_paired(bonds):
    '''
    :param bonds: in format {1: [[[2,1]], (2,1)], ...}
    :return: bonds in format {(c_1, c_2): [attr, length], ...}
    '''
    paired_bonds = {}
    for key, item in bonds.items():
        for i in item:
            paired_bonds.update({(key, i[0]): i[1::]})
    return paired_bonds


def bonds_to_one_way_dict(bonds):
    d = {}
    for i in bonds:
        if d.get(i[0]):
            d[i[0]].append(i[1::])
        else:
            d.update({i[0]: [i[1::]]})
    return d


def to_two_ways_bond2(one_way_bonds, with_attr=False):
    two_ways = {}
    if with_attr:
        for key, item in one_way_bonds.items():
            for i in sorted(item):
                if two_ways.get(key):
                    two_ways[key].append(i)
                else:
                    two_ways.update({key: [i]})
                if two_ways.get(i[0]):
                    two_ways[i[0]].append([key, i[1]])
                else:
                    two_ways.update({i[0]: [[key, i[1]]]})
        return two_ways
    for i in one_way_bonds:
        if two_ways.get(i[0]):
            two_ways.get(i[0]).append(i[1])
        else:
            two_ways.update({i[0]: [i[1]]})
        if two_ways.get(i[1]):
            two_ways.get(i[1]).append(i[0])
        else:
            two_ways.update({i[1]: [i[0]]})
    return two_ways


def mol2_to_notation(info_from_file, n_y=4, n_z=4, method='first', fr=None, sr=None):
    '''
    :param info_from_file: read_from_file tuple
    WARNING: may be keep Atoms without coordinates and Bonds with sections
    :return: chain in dictionary notation (see above)
    '''
    bonds, atoms = info_from_file
    positions_copy = copy.deepcopy(atoms)
    notation = {}
    bonds = bonds_to_one_way_dict(prepare_bonds(bonds))
    bonds2 = to_two_ways_bond2(bonds, with_attr=True)
    for key, item in atoms.items():
        cur_p = positions_copy.pop(key).position()
        connected = [i[0] for i in bonds2[key]]
        if n_y == None or n_z == None:
            basis = find_basis(cur_p, [atoms[i].position() for i in connected], method=method)
            # print(basis)
            notation.update({key: [list([[i, find_section(cur_p, atoms[i].position(), basis0=basis, method=method)]
                                         for i in connected]), basis]})
        else:
            basis = find_basis(cur_p, [atoms[i].position() for i in connected], n_y=n_y, n_z=n_z, method=method, fr=fr,sr=sr)
            # print(basis)
        atoms[key].set_orientation(basis)

        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position(), basis0=basis, n_y=n_y, n_z=n_z, method=method, fr=fr, sr=sr)]
                                     for i in connected]), basis]})
    for key, item in bonds.items():
        for i in range(len(item)):
            bonds[key][i].insert(1, np.linalg.norm(atoms[key].position()-atoms[item[i][0]].position()))
    return notation, bonds


def dimensional_structure(notation, n_y=None, n_z=None, fr=None, sr=None, method='first'):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    with length
    :return: xyz-info
    '''
    bonds_l, lengths = notation
    first_atom = min(bonds_l.keys())
    dim_structure = {first_atom: np.array([0, 0, 0])}#, np.array([0, 0])]}
    p = bonds_l[first_atom]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, first_atom)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds, basis = p.pop(0)
        for i in bonds:  # build bonds for cur_key atom
            if not (i[0] in dim_structure):  # if we don't have position:
                if method=='first':
                    if n_y!=None and n_z!=None:
                        coord = rotate_by_basis(pp[i[1]], basis[0], basis[1], n_y=n_y, n_z=n_z)*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                                                 else lengths[(i[0], cur_key)][0]) \
                                + dim_structure[cur_key]
                    else:
                        coord = rotate_by_basis(pp[i[1]], basis[0], basis[1])*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                           else lengths[(i[0], cur_key)][0])\
                            + dim_structure[cur_key]
                elif method == 'incline':
                    if n_y!=None and n_z!=None:
                        if fr!=None and sr!=None:
                            coord = rotate_non_perpendicular(pp[i[1]], basis[0], basis[1], n_y=n_y, n_z=n_z, fr=fr, sr=sr) * (
                                lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                else lengths[(i[0], cur_key)][0]) \
                                    + dim_structure[cur_key]
                        else:
                            coord = rotate_non_perpendicular(pp[i[1]], basis[0], basis[1], n_y=n_y, n_z=n_z)*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                                                 else lengths[(i[0], cur_key)][0]) \
                                + dim_structure[cur_key]
                    else:
                        coord = rotate_non_perpendicular(pp[i[1]], basis[0], basis[1])*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                           else lengths[(i[0], cur_key)][0])\
                            + dim_structure[cur_key]
                elif method == 'ten':
                    coord = rotate_ten_vars(pp[i[1]], basis)*(lengths[(cur_key, i[0])][0] if cur_key < i[0]
                                                                                                 else lengths[(i[0], cur_key)][0]) \
                                + dim_structure[cur_key]

                dim_structure.update({i[0]: coord})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                p.append(poper)
    return dim_structure


def write_mol2_file(file_name, atoms, positions, bonds):
    '''
    :param file_name:
    :param names: chemical elements names
    :param positions: xyz - coordinates
    :param bonds: one way bonds
    :return: None, void write function
    '''
    with open(file_name, 'w') as f1:
        f1.write('@<TRIPOS>MOLECULE\n')
        f1.write('some info\n')
        f1.write(str(len(atoms))+'\t'+str(len(bonds))+'\n\n')
        f1.write('@<TRIPOS>ATOM\n')
        for num, key in atoms.items():
            f1.write("\t{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(num, key.name, str(positions[num][0]),
                                                                            str(positions[num][1]), str(positions[num][2]),
                                                                            key.name_i, key.i1, key.i2, key.i3))
        f1.write('@<TRIPOS>BOND\n')

        for k, num in enumerate(bonds.items()):
            num, i = num
            # print(i, 't',attrs.get(tuple([i[0], i[1]])), attrs.get(tuple([i[1], i[0]])))
            f1.write("\t{0}\t{1}\t{2}\t{3}\n".format(str(k+1), str(num[0]), str(num[1]), str(i[1])))


##############################Mol2_preparations###########################################

def to_one_way_bonds(two_way_bonds):
    one_way = []
    for key, item in two_way_bonds.items():
        ix = 0
        if isinstance(item[0], (int, float)):
            ix = 1
        for i in item[ix]:
            if [i[0], key] not in one_way:
                one_way.append([key, i[0]])
    return one_way


def to_two_ways_bond(one_way_bonds, with_attr=False):
    two_ways = {}
    print(sorted(one_way_bonds))
    if with_attr:
        for i in one_way_bonds:
            if two_ways.get(i[0]):
                two_ways.get(i[0]).append((i[1], i[2]))
            else:
                two_ways.update({i[0]: [(i[1], i[2])]})
            if two_ways.get(i[1]):
                two_ways.get(i[1]).append((i[0], i[2]))
            else:
                two_ways.update({i[1]: [(i[0], i[2])]})
        return two_ways

    for i in one_way_bonds:
        if two_ways.get(i[0]):
            two_ways.get(i[0]).append(i[1])
        else:
            two_ways.update({i[0]: [i[1]]})
        if two_ways.get(i[1]):
            two_ways.get(i[1]).append(i[0])
        else:
            two_ways.update({i[1]: [i[0]]})
    return two_ways

def molecular_divider(atoms, bonds):
    '''
    :param atoms: bonds, atoms = xyz_names_bonds()
    :return: {atom_num: mol_num} (prob. later: div_atoms - dictionary like {mol_num: [atoms_of_this_molecule], ...},)
     [lig_as] - list of dictionaries with atoms of structure
     [lig_bs] - list of dictionaries with bonds of structure
    '''
    atom_names = set([j.i2 for _, j in atoms.items()])
    div_atoms, ligs_atoms, ligs_bonds = {}, {}, {}
    for ix, atom in enumerate(atom_names):
        lig_as = dict(filter(lambda x: x[1].i2 == atom, atoms.items()))
        lig_bs = list(filter(lambda x: x[0] in lig_as.keys() or x[1] in lig_as.keys(), bonds))
        for lig in lig_as:
            div_atoms.update({lig: ix})
        ligs_atoms.update({ix: lig_as})
        ligs_bonds.update({ix: lig_bs})
    return div_atoms, ligs_atoms, ligs_bonds


def get_notation_many_mols(atoms, bonds):
    '''
    :param atoms:
    :param bonds: bonds, atoms = xyz_names_bonds()
    :return: dictionary of notations for every molecule
    for atoms and for bonds separatedly
    '''
    div_atoms, ligs_as, lig_bonds = molecular_divider(atoms, bonds)
    not_atoms, not_bonds = {}, {}
    for i in ligs_as.keys():
        lig_as = ligs_as[i]
        lig_bs = lig_bonds[i]
        if lig_bs != []:
            ln = mol2_to_notation([lig_bs, lig_as])
            not_atoms.update({i: ln[0]})
            not_bonds.update({i: ln[1]})
        else:
            not_bonds.update({i: []})
            not_atoms.update({i: [ixx for ixx in lig_as.keys()][0]})
    return not_atoms, not_bonds


if __name__ == '__main__':
    '''Test: read from mol2-file by atoms_and_bonds() - function
    allows us get information about every atom;
    xyz_names_bonds()- function
    '''

    name = 'vacuum_cation_singlet_Fe_full'
    bs, ass = xyz_names_bonds(name + '.mol2')
    # atoms_info = atoms_and_bonds(name + '.mol2')
    # print(atoms_info)
    div_atoms, ligs_as, lig_bonds = molecular_divider(ass, bs)
    atoms_notation, bonds_notation = get_notation_many_mols(ass, bs)

    # print(div_atoms)
    atoms = set([j.i2 for _, j in ass.items()])
    nots = []
    div_atoms = {}
    for ix, atom in enumerate(atoms):
        lig_as = dict(filter(lambda x: x[1].i2 == atom, ass.items()))
        lig_bs = list(filter(lambda x: x[0] in lig_as.keys() or x[1] in lig_as.keys(), bs))
        if lig_bs != []:
            ln = mol2_to_notation([lig_bs, lig_as])
            nots.append(ln)
            # dd = dimensional_structure([ln[0], bonds_of_paired(ln[1])])
        else:
            nots.append([lig_as, []])
        for lig_one in lig_as:
            div_atoms.update({lig_one: ix})
        #     print(lig_as)
    # print([len(i[0]) for i in nots])
    finder = sorted(nots, key=lambda x: len(x[0]))
    mol_lengths = [len(i[0]) for i in finder]
    num_mols = len(finder)
    connecters = []
    for i in range(num_mols - 1):
        cur_atoms = set(finder[i][0].keys())
        meta_dists = []
        distances = []
        for k in cur_atoms:
            for key, ik in ass.items():
                distances.append((np.linalg.norm(ass[k].position()-ik.position()), key, k))
        connecteds = sorted(distances)
        fuller = {j: 0 for j in range(num_mols) if j != i}
        for j in connecteds[mol_lengths[i]::]:
            if div_atoms[j[1]] != i and div_atoms[j[1]] != div_atoms[j[2]] and fuller[div_atoms[j[1]]] < 2:
                connecters.append((j[2], j[1]))
                fuller[div_atoms[j[1]]] += 1
    print(connecters)
    zero_bonds = {}
    for ix, i in enumerate(connecters):
        c1, c2 = sorted(i)
        p1 = ass[i[0]].position()
        p2 = ass[i[1]].position()
        zero_bonds.update({ix: Bond(c1, c2, '0', length=np.linalg.norm(p1-p2), sections=[find_section(p1, p2), find_section(p2, p1)])})

                # print(ass[k])
        # print(finder[0][0].keys())
        # finder_c = list(filter(lambda x: set(x[0].keys()).intersection(cur_atoms) == {}, finder))
        # print(finder_c)

        # find_two_nearest_atoms
        # for j in cur_atoms.keys():
        #     finder_c.pop(j)
        # print(cur_atoms)



    # write_mol2_file("My_one_atom.mol2", lig_as, dd, bonds=bonds_of_paired(ln[1]))
    # (xyz_names_bonds(name + '.mol2'))
    # ln = mol2_to_notation(xyz_names_bonds(name + '.mol2'))
    # print(ln)
    # paired = bonds_of_paired(ln[1])
    # dim_structure = dimensional_structure([ln[0], paired])
    # write_mol2_file('My_'+name+'.mol2', atoms_info, dim_structure, bonds=paired)
