import copy
import numpy as np
from numpy import arctan2
from mol2_chain_q import bonds_of_paired



def cartesian_to_spherical(vector):
    r_xy = vector[0] ** 2 + vector[1] ** 2
    theta = arctan2(vector[1], vector[0])
    phi = arctan2(vector[2], r_xy ** 0.5)
    return theta, phi


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


def check(notation, dim_structure_reduced, eps=0.01):
    forces = 1
    forces_next = 0
    notation_l = notation.notation

    while abs(forces_next - forces) > eps**3:
        forces_next, forces = 0, forces_next
        lengths = bonds_of_paired(notation.bonds)
        for atom in dim_structure_reduced:
            force = np.array([0, 0, 0])
            for bond in notation_l[atom][0]:
                if bond[0] in dim_structure_reduced:
                    length = lengths.get(tuple(sorted([atom, bond[0]])))[0]
                    force = force + dim_structure_reduced[bond[0]]-notation.divider.scube[bond[1]]*length-dim_structure_reduced[atom]
            n_f = np.linalg.norm(force)
            # if n_f>0.05: print('n_f', n_f)
            forces_next += n_f
            dim_structure_reduced[atom] = dim_structure_reduced[atom] + eps*force
        # print(atom, forces_next)
    return dim_structure_reduced

def dimensional_structure(notation, relax=True):
    '''
    :param notation: Notation #TODO write more
    with length
    :return: xyz-info
    '''
    div = notation.divider.scube
    lengths = bonds_of_paired(notation.bonds)
    bonds_l = copy.deepcopy(notation.notation) # warning - it was error in other dim_structure builders
    print(bonds_l)
    first_atom = min(bonds_l.keys())
    dim_structure = {first_atom: np.array([0, 0, 0])}
    p = [list(bonds_l[first_atom].keys())]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, first_atom)  # p[0] - current atom, p[1] - bonds
    p = [p]
    print(p[0])
    while len(p) != 0:
        cur_key, bonds = p.pop(0)
        print(cur_key)
        print(bonds)
        for i in bonds:  # build bonds for cur_key atom
            print(i)
            if not (i[0] in dim_structure):  # if we don't have position:
                coord = div[i[1]]*(lengths.get(tuple(sorted([cur_key, i[0]])))[0]) + dim_structure[cur_key]
                dim_structure.update({i[0]: coord})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                ix = -1
                while ix < len(poper[1])-1:
                    ix += 1
                    if poper[1][ix][0] == cur_key:
                        poper[1].pop(ix)
                        break
                p.append(poper)
            else:
                if relax: dim_structure = check(notation, dim_structure)
    return dim_structure

def relaxing(notation, lengths, dim_structure):
    scube = notation.divider.scube
    for i, j in notation.items():
        for k in j[0]:#i, k[0] elements considered
            delta_length = np.linalg.norm(dim_structure[k[0]]-(dim_structure[i]+scube[k[1]]*(lengths.get((i, k[0]), lengths.get(k[0], i))[0])))
            if i > k[0] and (delta_length) > 0.09:# and i < k[0]:
                print(i, k[0], round(delta_length, 3), lengths.get((i, k[0]), lengths.get(k[0], i))[1])
    return 0



if __name__ == '__main__':
    # name_sh = '4c-Mn-OMe'
    # names = ['Caffein', 'Naphthalene', 'Phenol', '4c-Mn-OMe', '3-MnH2', '2-Mn-OtBu', 'Mn-deprot-Mn-bare', 'Heroin_2']
    names = ['4c-Mn-OMe', '3-MnH2', '2-Mn-OtBu', 'Mn-deprot-Mn-bare']
    # names = ['H2']
    for name_sh in names:
        file_name = './tmp/'+name_sh+'_opt'
        name = name_sh+'_opt'
        # file_name = name_sh
        # file_name = name
        n_param = 6

        # bs, ass = xyz_names_bonds(name + '.mol2')


        # atoms_info = atoms_and_bonds(file_name + '.mol2')
        # ln = Notation(n=n_param, info_from_file=xyz_names_bonds(file_name + '.mol2'))
        # paired = bonds_of_paired(ln.bonds)
        # dim_structure = dimensional_structure(ln, relax=True)
        # write_mol2_file('My_' + name + '_' + 'q0.mol2', atoms_info, dim_structure, bonds=paired)
        #
        # file_name_new = 'My_' + name + '_' + 'q0'
        # atoms_info2 = atoms_and_bonds(file_name_new + '.mol2')
        # ln2 = Notation(n=n_param, info_from_file=xyz_names_bonds(file_name_new + '.mol2'))
        #
        # flags = []
        # for i in ln.notation.keys():
        #     flags.append(ln.notation[i].sort()==ln2.notation[i].sort())
        #     if ln.notation[i].sort() != ln2.notation[i].sort():
        #         print(ln.notation[i], ln2.notation[i])
        # print(name_sh, len(flags) - sum(flags), 'errors')
        #
        # flags = []
        # for i in ln.bonds.keys():
        #     b1, b2 = sorted(ln.bonds[i]), sorted(ln2.bonds[i])
        #     for j in range(len(b1)):
        #         flags.append(b1[j][1] == b2[j][1])
        #         if b1[j][1] != b2[j][1]:
        #             print(b1[j][1], b2[j][1])
        # print(name_sh, len(flags) - sum(flags), 'length errors')

        # print(ln.notation == ln2.notation)

        # dim_structure = dimensional_structure(ln2, relax=True)
        # write_mol2_file('My_' + name + '_' + 'q1.mol2', atoms_info, dim_structure, bonds=paired)
