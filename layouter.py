import copy
import numpy as np
from mol2_chain_q import atoms_and_bonds,  bonds_of_paired, mol2_to_notation, dimensional_structure, xyz_names_bonds, write_mol2_file
from quadro_with_rotate import scube, find_section, anti_scube
from numpy import arctan2, pi

def cartesian_to_spherical(vector):
    r_xy = vector[0] ** 2 + vector[1] ** 2
    theta = arctan2(vector[1], vector[0])
    phi = arctan2(vector[2], r_xy ** 0.5)
    return theta, phi

def dimensional_structure(notation, **kwargs):
    '''
    :param notation: Notation with first atom with unique basis for every bond
    with length
    :return: xyz-info
    '''
    bonds_l, lengths = copy.deepcopy(notation) #warning - it was error in other dim_structure builders
    first_atom = min(bonds_l.keys())
    dim_structure = {first_atom: np.array([0, 0, 0])}#, np.array([0, 0])]}
    p = bonds_l[first_atom]
    bonds_copy = copy.deepcopy(bonds_l)
    p.insert(0, first_atom)  # p[0] - current atom, p[1] - bonds, p[2] - basis of p[0] atom
    p = [p]
    while len(p) != 0:
        cur_key, bonds = p.pop(0)
        # print(cur_key, bonds)
        for i in bonds:  # build bonds f    or cur_key atom
            if not (i[0] in dim_structure):  # if we don't have position:
                # print((lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]))
                coord = scube[i[1]]*(lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]) + dim_structure[cur_key]
                dim_structure.update({i[0]: coord})
                poper = bonds_copy.pop(i[0])
                poper.insert(0, i[0])
                # poper.append(i[0])
                ix = -1
                while ix < len(poper[1])-1:
                    ix += 1
                    if poper[1][ix][0] == cur_key:
                        poper[1].pop(ix)
                        break
                p.append(poper)
            else:
                print('cycle:', cur_key, i[0])
                # print((lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]))

                # coord = scube[i[1]] * (lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]) + \
                #         dim_structure[cur_key]
                # delta_norm = np.linalg.norm(np.array(dim_structure[i[0]])-np.array(coord))
                # delta_phi, delta_psi = cartesian_to_spherical(scube[find_section(dim_structure[i[0]], dim_structure[cur_key])]-scube[i[1]])
                #
                # print(delta_norm)
                #
                # print(delta_phi*180/pi, delta_psi*180/pi)
                # forces = []
                # print(dim_structure.keys())
                # for i in dim_structure.keys():
                    # for j in notation[0][i]:
                    #     print('j', j)

                # for i in notation[0].keys():
                #     for _, j in notation[0].items():
                #         print(j)
                #     print(i)

                # print(cartesian_to_spherical(scube[find_section(i[0], dim_structure[cur_key])]),
                #       cartesian_to_spherical(scube[i[1]]))
                # print()

                # if np.linalg.norm(np.array(dim_structure[i[0]])-np.array(coord)) > 0.1:

                #     print(np.linalg.norm(dim_structure[i]))
                # print(np.linalg.norm(dim_structure[i[0]]-dim_structure[cur_key]))
                # print(dim_structure[cur_key])
                # print(dim_structure[i[0]])
    return dim_structure

def relaxing(notation, lengths, dim_structure):
    print(notation)
    print(dim_structure)
    forces = {}
    for i, j in notation.items():
        for k in j[0]:#i, k[0] elements considered
            delta_length = np.linalg.norm(dim_structure[k[0]]-(dim_structure[i]+scube[k[1]]*(lengths.get(tuple([i, k[0]]), lengths.get(tuple([k[0], i])))[0])))
            delta_psi, delta_phi = 0, 0 #TODO
            if (delta_length) > 0.1:# and i < k[0]:
                print(i, k[0], delta_length, lengths.get(tuple([i, k[0]]), lengths.get(tuple([k[0], i])))[1])

    return 0


if __name__ == '__main__':
    file_name = './tmp/Phenol_opt'
    name = 'Phenol_opt'
    print(anti_scube())

    # bs, ass = xyz_names_bonds(name + '.mol2')

    atoms_info = atoms_and_bonds(file_name + '.mol2')
    ln = mol2_to_notation(xyz_names_bonds(file_name + '.mol2'))
    print(ln)
    paired = bonds_of_paired(ln[1])
    dim_structure = dimensional_structure([ln[0], paired])
    relaxing(ln[0], paired, dim_structure)
    # print(dim_structure)
    write_mol2_file('My_' + name + '_' + 'q.mol2', atoms_info, dim_structure, bonds=paired)