from Ordered_edge import Ordered_edge
from Ordered_node import Ordered_node
class Ordered_Graph():
    def __init__(self, edges_list, dict_list):
        self.stack_nodes = {}
        self.set_of_nodes_nums = set()
        for i in edges_list:
            if not i[0] in self.set_of_nodes_nums:
                self.set_of_nodes_nums.add(i[0])
                self.stack_nodes[i[0]] = Ordered_node(dict_list[i[0]])
            if not i[1] in self.set_of_nodes_nums:
                self.set_of_nodes_nums.add(i[1])
                self.stack_nodes[i[1]] = Ordered_node(dict_list[i[1]])
            edge = Ordered_edge(self.stack_nodes[i[0]], self.stack_nodes[i[1]])
            self.stack_nodes[i[0]].order_edges.append(edge)
            self.stack_nodes[i[1]].order_edges.append(edge)
            # self.stack_nodes[i[0]].order_edges.append(self.stack_nodes[i[1]])
            # self.stack_nodes[i[1]].order_edges.append(self.stack_nodes[i[0]])

        # for item, key in dict_list:
        #     node = Ordered_node()

if __name__ == '__main__':
    gr = Ordered_Graph(None, [[1, 2],[1,3]], {1: 'C', 2:'H', 3:'H'})
    # print(gr.stack_nodes)
    # for i in gr.stack_nodes.items():
    #     print(i[1].element)
