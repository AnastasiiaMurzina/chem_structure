import matplotlib.pyplot as plt
import numpy as np
import cmath
from cmath import pi

def get_main_orient_part(structure):
    i = 0
    while structure[i][0] != structure[i + 1][0]:
        i += 1
    return structure[i][0]

def list_of_many_contrains_atoms(structure):
    ss = set()
    for i in structure:
        ss.add(i[0])
        ss.add(i[1])
    print(ss)

if __name__ == "__main__":
    '''Description of structure: list_of_bonds'''
    test_structure = [[1, 2, 1], [1, 3, 1], [1, 4, 1], [4, 5, 1], [4, 6, 1]] # last num is count of bonds in one
    t_s = {"bonds": tuple(test_structure),
           "order": {1: [1, 2], 4: [2, 1]} # [2, 1] rotate 1-4 pi/count_v
	    }
    #for i in range(len([i[0] if i[0]==1 for i in test_structure])):
	
    dict_structure = {1: "C", 2: "A1_group", 3: "A2_group", 4: "C", 5: "A3_group", 6: "A4_group"}
    max_possible_atom_bonds = {"C": 4, "A1_group": 1, "A2_group": 1,
                               "A3_group": 1, "A4_group": 1,
                               "A5_group": 1}
    distance_dict = {("C", "A1_group"): 1,
                     ("C", "A2_group"): 2,
                     ("C", "A3_group"): 2,
                     ("C", "A4_group"): 1,
                     ("C", "C"): 3}
    energy_bond = {("C", "A1_group"): 10,
                    ("C", "A2_group"): 15,
                     ("C", "A3_group"): 20,
                     ("C", "A4_group"): 10}
    structure_with_angles = {}

    cur_point = [0, 0]
    for i in t_s["order"]:
        pass
    # print(t_s["order"][1])
    # list_of_many_contrains_atoms(test_structure)
    # test_structure.append([4,7])
    # dict_structure.update({7: "A1_group"})

    coord_first = [0, 0]
    plt.scatter(coord_first[0], coord_first[1], marker='o')
    step_bounds = []
    for i in test_structure: # check bonds and save them
        if i[0] == 1:
            step_bounds.append(i[1])
        elif i[1] == 1:
            step_bounds.append(i[0])

    for index, i in enumerate(step_bounds): # set angles to bonds
        angle = index*2*pi/max_possible_atom_bonds[dict_structure[1]]
        endx = coord_first[0] + distance_dict[dict_structure[1], dict_structure[i]]*cmath.cos(angle)
        endy = coord_first[1] + distance_dict[dict_structure[1], dict_structure[i]]*cmath.sin(angle)
        structure_with_angles.update({(1, i): angle}) # * 180/(pi) to degrees
        plt.plot([coord_first[0], endx], [coord_first[1], endy], 'r-')
        # print(dict_structure[i])
        plt.annotate(
                    dict_structure[i],
                    xy=((endx.real), (endy.real)), xytext=(-20, 20),
                    textcoords='offset points', ha='right', va='bottom',
                    bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
                    arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    new_c = [endx, endy]
    step_bounds = []
    print(structure_with_angles)

    for i in test_structure:
        if i[0] == 4:
            step_bounds.append(i[1])
        elif i[1] == 4:
            angle = structure_with_angles[i[0], 4]
            step_bounds.append(i[0])

    for index, i in enumerate(step_bounds):
        angle2 = (angle + pi) + index * 2 * pi / max_possible_atom_bonds[dict_structure[4]]
        print(index, i, dict_structure[i], angle2)
        # angle2*=pi/180
        endx = new_c[0] + distance_dict[dict_structure[4], dict_structure[i]] * cmath.cos(angle2)
        endy = new_c[1] + distance_dict[dict_structure[4], dict_structure[i]] * cmath.sin(angle2)
        if not((i,4) in structure_with_angles.keys()): structure_with_angles.update({(4,i): angle2})
        plt.plot([new_c[0], endx], [new_c[1], endy], 'b-')
        plt.annotate(
            dict_structure[i],
            xy=((endx.real), (endy.real)), xytext=(-20, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'))
    plt.scatter(new_c[0], new_c[1], marker='o')
    print(structure_with_angles)

plt.show()

