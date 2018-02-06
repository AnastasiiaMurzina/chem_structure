import networkx as nx
import matplotlib.pyplot as plt
import random
dict_structure = {1: "C", 2: "A1_group", 3: "A2_group", 4: "C", 5: "A3_group", 6: "A4_group",
                  7: "B1_group", 8: "B2_group", 9: "B3_group", 10: "B4_group",11: "B5_group",
                  12: "C1_group", 13: "C2_group"}
def penten_structure():
    penten = [(1,4), (2,4), (3,4), (4,5), (5,6), (5,7), (7,15), (7,8), (8,9), (8,10), (8,11),(11,12),(11,13),(11,14)]
    dict_penten = {1: "H", 2: "H", 3: "H", 4: "C", 5:"C", 6: "H", 7: "C", 8: "C", 9: "H", 10: "H", 11:"C", 12: "H", 13: "H", 14: "H",15: "H"}
    Penten = nx.Graph()
    Penten.add_edges_from(penten)
    for i in Penten.edges(data=True):
        i[2]['multiple_bond'] = 1
        i[2]['part_of_cycle'] = False
    Penten.get_edge_data(5,7)['multiple_bond'] = 2
    return Penten, dict_penten

def show_graph(graph):
    plt.subplot(121)
    nx.draw(graph, with_labels=True, font_weight='bold')
    # plt.subplot(122)
    # nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')
    plt.show()

def cyclooctan_structure():
    cyclooctan = []
    for i in range(8):
        cyclooctan.append((i+1, max((i+2) % 9, 1)))
    for i in range(6):
        cyclooctan.append((i+2, i+9))
        cyclooctan.append((i+2, 20-i))
    cyclooctan.append((1,21))
    cyclooctan.append((8,22))
    dict_cyclooctan = {}
    for i in range(8):
        dict_cyclooctan[i + 1] = 'C'
    for i in range(14):
        dict_cyclooctan[i + 9] = 'H'
    Cyclooctan = nx.Graph()
    Cyclooctan.add_edges_from(cyclooctan)
    Cyclooctan.get_edge_data(1,8)['multiple_bond'] = 2
    for i in range(len(cyclooctan)):
        for j in Cyclooctan.edges(i+1, data=True):
            j[2]['part_of_cycle'] = False
        if len(Cyclooctan.edges(i+1)) == 4 or len(Cyclooctan.edges(i+1)) == 1:
            for j in Cyclooctan.edges(i+1, data=True):
                j[2]['multiple_bond'] = 1
    for i in range(8):
        Cyclooctan.get_edge_data(i+1, max((i+2) % 9, 1))['part_of_cycle'] = True
    # show_graph(Cyclooctan)
    return Cyclooctan, dict_cyclooctan
# cyclooctan_structure()
show_graph(penten_structure()[0])
def bds(structure, num_zeros_of_go, first_node, current_step = [0, 0]):
    G, d = structure
    if current_step == [0, 0]:
        stek = [first_node]
        visited = [num_zeros_of_go]
    else:
        stek, visited = current_step
    stack = []
    for i in stek:
        for j in G.edges(i, data=True):
            if not j[1] in visited:
                # print(j)
                # print(d[j[1]])
                stack.append(j[1])
        visited.append(i)
    return stack, stek# visited
# bds(penten_structure(), 1, first_node=4)
# bds(penten_structure(), 1, 4, bds(penten_structure(), 1, first_node=4))

def get_graph():
    G = nx.Graph()
    G.add_node(1)
    G.add_edge(1, 2, multiple_bond=1,
               part_of_cycle=random.choice([True, False]),
               count_of_combination_properties=0)
    for i in range(15):
        G.add_edge(i + 2, i + 3, multiple_bond=random.randint(1, 2),
                   part_of_cycle=random.choice([True, False]),
                   count_of_combination_properties=0)
    return G

def get_alone_sub_around_cycle(graph, i):
    array_alones = graph.edges(i)

def set_alkene_isomerism(graph):
    step = [0, 0]
    # G, d = penten_structure()
    G, d = cyclooctan_structure()
    G = nx.Graph(G)
    #alkene isomers of double bound without cycles at all
    for i in G.edges(data=True):
        if i[2]['multiple_bond'] == 2:
            if not i[2]['part_of_cycle']:
                i[2]['alkene_cycle'] = False
                # print(i)
                # print(bds([G, d], i[0], i[1], current_step = [0, 0]))
                # print(bds([G, d], i[1], i[0], current_step = [0, 0]))
                bds1 = bds([G,d], i[0], i[1], current_step = [0, 0])
                bds2 = bds([G,d], i[1], i[0], current_step = [0, 0])

                compare1 = [d[j] for j in bds1[0]]
                compare2 = [d[j] for j in bds2[0]]
                # very IMPORTANT do recursion
                if compare1[0]!=compare1[1] and compare2[0]!=compare2[1]:
                    i[2]['alkene_noncycle'] = True
                    print('it is alkene noncycle isomerism', i)
                else:
                    i[2]['alkene_noncycle'] = False
            else:
                i[2]['alkene_noncycle'] = False
                lonely_elem1 = [j for j in G.edges(i[0], data=True) if not j[2]['part_of_cycle']]
                lonely_elem2 = [j for j in G.edges(i[1], data=True) if not j[2]['part_of_cycle']]
# BIG question here: only hydrogen?????
                if d[lonely_elem1[0][1]] == d[lonely_elem2[0][1]] == 'H':
                    i[2]['alkene_cycle'] = True
                    print('it is alkene cycle isomerism', i)
        else:
            i[2]['alkene_noncycle'] = False
            i[2]['alkene_cycle'] = False
    # for i in G.edges(data=True):
    #     print(i)

# set_alkene_isomerism(None)
G = get_graph()
for i,j,data in G.edges(data=True):
    # print(data)
    if data['multiple_bond'] == 2:
        if data['part_of_cycle']:
            element1_cis_trans_cyclo = [x for x in G.edges(i, data=True) if x[2]["part_of_cycle"] == False]
            element2_cis_trans_cyclo = [x for x in G.edges(j, data=True) if x[2]["part_of_cycle"] == False]
            # print("CIS", element1_cis_trans_cyclo, element2_cis_trans_cyclo)
            # if element1_cis_trans_cyclo != [] and dict_structure[element1_cis_trans_cyclo[0]] == "A4_group" and element2_cis_trans_cyclo != [] and dict_structure[element2_cis_trans_cyclo[0]] == "A5_group":
            #     print("alkene_cyclo_isomer")


        # print(G.edges(i, data = True))
        # print(G.edges(j))
        # print('alkene isomer')
    # print(n, G[n[0]][n[1]]['multiple_bond'])
# e=(2,3)

# print(random.randint(0,2))
# G.add_edge(*e)
# print(G.get_edge_data(2,3))
# print(G[1][2]['multiple_bond'])
# G.add_edge(2,3)
# plt.subplot(121)
# nx.draw(Penten, with_labels=True, font_weight='bold')
# plt.subplot(122)
# nx.draw_shell(G, nlist=[range(5, 10), range(5)], with_labels=True, font_weight='bold')
# plt.show()
