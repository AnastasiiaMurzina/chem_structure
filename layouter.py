import copy
import numpy as np
from numpy import arctan2


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
            f1.write("\t{0}\t{1}\t{2}\t{3}\n".format(str(k+1), str(num[0]), str(num[1]), str(i[1])))


def check(notation, dim_structure_reduced, eps=0.01):
    forces = 1
    forces_next = 0
    while abs(forces_next - forces) > eps: #**3
        forces_next, forces = 0, forces_next
        for atom in dim_structure_reduced:
            force = np.array([0, 0, 0])
            for bond in notation.bonds[atom].keys():
                if bond in dim_structure_reduced:
                    l = notation.bonds[atom][bond].length
                    s = notation.bonds[atom][bond].section
                    d = dim_structure_reduced[atom]+notation.divider.scube[s]*l-dim_structure_reduced[bond]
                    force = force + d
            n_f = np.linalg.norm(force)
            forces_next += n_f
            dim_structure_reduced[atom] = dim_structure_reduced[atom] - eps*force
    return dim_structure_reduced

def dimensional_structure(notation, relax=True, eps=0.01):
    '''
    :param notation: Notation bfs
    :return: xyz dictionary
    '''
    div = notation.divider.scube
    bonds_l = copy.deepcopy(notation.notation) # warning - it was error in other dim_structure builders
    first_atom = min(bonds_l.keys())
    dim_structure = {first_atom: np.array([0, 0, 0])}
    p = [[first_atom, list(bonds_l[first_atom].keys())]] # p[0] - current atom, p[1] - bonds
    while len(p) != 0:
        cur_key, bonds = p.pop(0)
        for i in bonds:  # build bonds for cur_key atom
            if not (i in dim_structure):  # if we don't have position:
                s = notation.notation[cur_key][i]
                coord = div[s]*notation.bonds[cur_key][i].length + dim_structure[cur_key]
                dim_structure.update({i: coord})
                p.append([i, list(bonds_l.pop(i).keys())])
                if relax: dim_structure = check(notation, dim_structure, eps=eps)
    return dim_structure


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
