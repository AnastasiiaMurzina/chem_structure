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

def all_neighbor_atoms(structure, num):
    neighbors = []
    bound = max_possible_atom_bonds[dict_structure[num]]
    for i in structure:
        if i[0] == num:
            neighbors.append(i[1])
            bound -= i[2]
        if i[1] == num:
            neighbors.append(i[0])
            bound -= i[2]
        if bound == 0:
            break
    return neighbors

def add_alkene_isomerism(structure, dict_structure):
    for i in structure:
        if i[2] == 2 and max_possible_atom_bonds[dict_structure[i[1]]]>3 and max_possible_atom_bonds[dict_structure[i[0]]]>3:
            this_node, next_node, count = i[0], i[1], i[2]
            this_node_first, this_node_second = [j for j in all_neighbor_atoms(structure, this_node) if j!= next_node]
            next_node_first, next_node_second = [j for j in all_neighbor_atoms(structure, next_node) if j!= this_node]
            # or there is a different on the next steps!!!!!
            if this_node_first!=this_node_second or next_node_first!=next_node_second:
                print("this bond", this_node, "-", next_node, " may be in two statements: cis- and tranc- (0 or 1)")
        # also possible in  cyclic part of structure!!!

if __name__ == "__main__":
    '''Description of structure: list_of_bonds'''
    test_structure = [[1, 2, 1], [1, 3, 1], [1, 4, 2], [4, 5, 1], [4, 6, 1]] # last num is count of bonds in one
    # t_s = {"bonds": tuple(test_structure),
    #        "order": {1: [1, 2], 4: [2, 1]}  # [2, 1] rotate 1-4 pi/count_v
	 #    }
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
    # structure_with_angles = {}
    add_alkene_isomerism(test_structure, dict_structure)



