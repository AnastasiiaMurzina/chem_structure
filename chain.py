from dividers.penta import get_penta_points, get_dictionary_coordinates, section_num_by_coords, show_points
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
count_bonds = {'H': 1, 'C': 4, 'O': 2}
eps_length = 0.001

def chain1():
    return [[1,2,1], [1,3,2], [1,4,6]], {1: 'C', 2: 'H', 3: 'H', 4: 'H'}

def chain2():
    return [[1,2,1], [2,3,5], [3,4,7]], {1: 'H', 2: 'O', 3: 'O', 4: 'H'}

def chain3():
    return [[1,2,1], [2,3,5], [3,4,7], [3,5, 10], [5,6,12], [3,7,3]],\
           {1: 'H', 2: 'O', 3: 'C', 4: 'H', 5: 'O', 6: 'H', 7: 'H'}

def show_with_bonds(bonds, dpoints, annotate=False, dictionary=None):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in bonds:
        ax.plot([dpoints[i[0]][0], dpoints[i[1]][0]],
                [dpoints[i[0]][1], dpoints[i[1]][1]],
                [dpoints[i[0]][2], dpoints[i[1]][2]])
        # print((dpoints[i[0]][0] - dpoints[i[1]][0]) ** 2 + (dpoints[i[0]][1] - dpoints[i[1]][1]) ** 2 +
        #       (dpoints[i[0]][2] - dpoints[i[1]][2]) ** 2)
    if annotate:
        for key, item in dpoints.items():
            ax.text(item[0]+0.05, item[1]+0.05, item[2]+0.05, dictionary[key])
    plt.show()

def show_points(dpoints, dict_structure):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for key, i in dpoints.items():
        ax.scatter(i[0], i[1], i[2])
        ax.text(i[0]+0.15, i[1]+0.15, i[2]+0.15, dict_structure[key])
    plt.show()


def dimetional_structure_from_list(list_bonds, dict_structure):

    first_node = [0, 0, 0]
    penta_points = get_penta_points(first_node)
    init_dict = get_dictionary_coordinates(penta_points)

    # list_structure, dict_structure = chain3()
    list_structure_copy = list_bonds.copy()
    queue = [1,]
    set_viewed_elements = set([])
    dict_structure_position = {1: np.array([0, 0, 0])}
    while len(queue) != 0:
        current_node = queue.pop()
        pops = []
        for num, i in enumerate(list_structure_copy):
            if i[0] == current_node and not i[0] in set_viewed_elements:
                dict_structure_position.update({i[1]: dict_structure_position[i[0]]+init_dict[i[2]]})
                queue.insert(0, i[1])
                pops.append(num)
        for i in pops.__reversed__():
            list_structure_copy.pop(i)
    return dict_structure_position, dict_structure

def write_to_file_coord_point(file, elements):
    with open(file, 'w') as f:
        for key, item in elements.items():
            position = dict_structure_position[key]
            string_arr = [item, position[0], position[1], position[2]]
            string = ' '.join([str(i) for i in string_arr])
            f.write(string+'\n')

def dictionary_for_checked_bond_count(dict_structure):
    dict_structure_copy = dict_structure.copy()
    for key, item in dict_structure_copy.items():
        dict_structure_copy[key] = count_bonds[item]
    return dict_structure_copy

def from_coordinates_to_list_bond(dict_structure_position, dict_structure):
    dict_structure_copy = dictionary_for_checked_bond_count(dict_structure)
    dict_structure_position_copy = dict_structure_position.copy()
    bonds_list = []
    for key, item in dict_structure_position.items():
        current_node = dict_structure_position_copy.pop(key)
        for i in range(dict_structure_copy[key]):
            for key2, item2 in dict_structure_position_copy.items():
                if abs(np.linalg.norm(current_node - item2) - 1) < eps_length and dict_structure_copy[key2] != 0:
                    bonds_list.append([key, key2, section_num_by_coords(current_node, item2)])
                    dict_structure_copy[key] -= 1
                    dict_structure_copy[key2] -= 1
                if dict_structure_copy[key] == 0: break
    return bonds_list

if __name__ == '__main__':
    struct = chain3()
    dict_structure_position, dict_structure = dimetional_structure_from_list(struct[0], struct[1])
    print(dict_structure_position)
    # show_with_bonds()
    print(from_coordinates_to_list_bond(dict_structure_position, dict_structure))

    # write_to_file_coord_point('output.txt', dict_structure, dict_structure_position)

    show_with_bonds(struct[0], dict_structure_position, annotate=True, dictionary=dict_structure)
    # show_with_bonds(list_structure, dict_structure_position)
    # show_points(dict_structure_position, dict_structure)
