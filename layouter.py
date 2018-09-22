import copy
import numpy as np
from mol2_chain_q import atoms_and_bonds,  bonds_of_paired, mol2_to_notation, dimensional_structure, xyz_names_bonds, write_mol2_file
from quadro_with_rotate import scube



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
        # print(cur_key, bonds)
        for i in bonds:  # build bonds f    or cur_key atom
            if not (i[0] in dim_structure):  # if we don't have position:
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
                coord = scube[i[1]] * (lengths.get(tuple([cur_key, i[0]]), lengths.get(tuple([i[0], cur_key])))[0]) + \
                        dim_structure[cur_key]
                print(np.linalg.norm(np.array(dim_structure[i[0]])-np.array(coord)))
                # print()

                # if np.linalg.norm(np.array(dim_structure[i[0]])-np.array(coord)) > 0.1:

                #     print(np.linalg.norm(dim_structure[i]))
                # print(np.linalg.norm(dim_structure[i[0]]-dim_structure[cur_key]))
                # print(dim_structure[cur_key])
                # print(dim_structure[i[0]])
    return dim_structure


if __name__ == '__main__':
    file_name = './tmp/Phenol_opt'
    name = 'Phenol_opt'

    # bs, ass = xyz_names_bonds(name + '.mol2')

    atoms_info = atoms_and_bonds(file_name + '.mol2')
    ln = mol2_to_notation(xyz_names_bonds(file_name + '.mol2'))
    print(ln)
    paired = bonds_of_paired(ln[1])
    dim_structure = dimensional_structure([ln[0], paired])
    # print(dim_structure)
    write_mol2_file('My_' + name + '_' + 'q.mol2', atoms_info, dim_structure, bonds=paired)