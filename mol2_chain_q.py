import numpy as np
import copy
from quadro_with_rotate import spherical_cube, find_section
from mol2_worker import xyz_names, xyz_names_bonds, Atom, atoms_and_bonds, Bond
from mopac_worker import get_energy_of_mol2
# from many_mols import molecular_divider, get_notation_many_mols, insert_zero_bonds

eps_length = 0.001
scube = spherical_cube(3)

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


def mol2_to_notation(info_from_file, **kwargs):
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
        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position())]
                                         for i in connected])]})
        notation.update({key: [list([[i, find_section(cur_p, atoms[i].position())]
                                     for i in connected])]})
    for key, item in bonds.items():
        for i in range(len(item)):
            bonds[key][i].insert(1, np.linalg.norm(atoms[key].position()-atoms[item[i][0]].position()))
    return notation, bonds


def dimensional_structure(notation, **kwargs):
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
        cur_key, bonds = p.pop(0)
        for i in bonds:  # build bonds for cur_key atom
            # print(i)
            # return
            if not (i[0] in dim_structure):  # if we don't have position:
                coord = scube[i[1]]*(lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]) + dim_structure[cur_key]

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



if __name__ == '__main__':
    '''Test: read from mol2-file by atoms_and_bonds() - function
    allows us get information about every atom;
    xyz_names_bonds()- function
    '''

    name = 'Cocaine'
    # bs, ass = xyz_names_bonds(name + '.mol2')
    atoms_info = atoms_and_bonds(name + '.mol2')
    ln = mol2_to_notation(xyz_names_bonds(name + '.mol2'))
    print(ln)
    paired = bonds_of_paired(ln[1])
    dim_structure = dimensional_structure([ln[0], paired])
    # print(dim_structure)
    write_mol2_file('My_'+name+'_'+'q.mol2', atoms_info, dim_structure, bonds=paired)

    ###################################Check_3_9######################################333
    import os, glob, subprocess
    def create_mol2s():
        g_dir = os.path.join(os.getcwd(), 'mols_dir', 'Molecules')
        r_dir = os.path.join(os.getcwd(), 'tmp', 'ico_first')
        for i in [x[0] for x in os.walk(g_dir)][1::]:
            cur_d = os.path.join(g_dir, i.split('/')[-1])
            cur_r_d = os.path.join(r_dir, os.path.basename(os.path.normpath(cur_d)))
            os.mkdir(cur_r_d)
            files = glob.glob(i + '/*_opt.xyz')
            for xyz_opt_file in files:
                # print(xyz_opt_file)
                mol2_file = os.path.join(r_dir, os.path.basename(os.path.normpath(xyz_opt_file)))
                subprocess.call(['babel', '-ixyz', xyz_opt_file, '-omol2', mol2_file[:-4:]+'.mol2'])

    def hist_erros(mth):
        r_dir = os.path.join(os.getcwd(), 'tmp', 'ico_first')
        mth = 'first'
        n_y, n_z = 4, 4
        # print(r_dir)
        graphic = []
        for name in glob.glob(r_dir + '/*_opt.mol2'):
            # print(name)
            exactly_name = os.path.basename(os.path.normpath(name))
            bs, ass = xyz_names_bonds(name)
            atoms_info = atoms_and_bonds(name)

            ln = mol2_to_notation(xyz_names_bonds(name), method=mth, n_y=n_y, n_z=n_z, fr=9, sr=13,
                                  r=0.6)  # , kwargs={'n_y': 5, 'n_z': 7})
            # print(ln)
            paired = bonds_of_paired(ln[1])
            dim_structure = dimensional_structure([ln[0], paired], method=mth, n_y=n_y, n_z=n_z, fr=9, sr=13,
                                                  r=0.6)  # ,kwargs={'n_y': 5, 'n_z': 7})
            dec_file = os.path.join(r_dir, 'My_' + exactly_name[:-5:] + '_' + mth + '.mol2')
            write_mol2_file(dec_file, atoms_info, dim_structure, bonds=paired)
            print(exactly_name)
            en_norm = get_energy_of_mol2(name)
            en_no_norm = get_energy_of_mol2(dec_file)
            print(en_no_norm, en_norm, 'cycle' if open(name, 'r').read().count('ar') >= 3 else '')
            graphic.append(abs((en_no_norm-en_norm)/en_norm))

        import matplotlib.pyplot as plt
        plt.hist(graphic)
        plt.xlabel('\delta energy / energy_of_original')
        plt.show()
        from mopac_worker import get_energy
        # print(exactly_name)

        # files = glob.glob(i + '/*.arc')

    # bs, ass = xyz_names_bonds(name + '.mol2')
    # atoms_info = atoms_and_bonds(name + '.mol2')
    # print(atoms_info)


    # write_mol2_file("My_one_atom.mol2", lig_as, dd, bonds=bonds_of_paired(ln[1]))
    # (xyz_names_bonds(name + '.mol2'))
    # ln = mol2_to_notation(xyz_names_bonds(name + '.mol2'), method=mth, n_y=n_y, n_z=n_z, fr=9, sr=13,
    #                       r=0.6)  # , kwargs={'n_y': 5, 'n_z': 7})
    # print(ln)
    # paired = bonds_of_paired(ln[1])
    # dim_structure = dimensional_structure([ln[0], paired], method=mth, n_y=n_y, n_z=n_z, fr=9, sr=13,
    #                                       r=0.6)  # ,kwargs={'n_y': 5, 'n_z': 7})
    # write_mol2_file('My_' + name + '_' + mth + '.mol2', atoms_info, dim_structure, bonds=paired)