from Ordered_node import Ordered_node
class Ordered_edge():

    rotation = 0
    length = 0
    energy = 0

    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2

    def set_rotation(self, degrees): # set zone by divider12
        self.rotation = degrees


    # def __del__(self):
    #     del self.node1
    #     del self.node2

    def node1_atoms(self):
        return self.node1.order_edges

    def node2_atoms(self):
        return self.node2.order_edges

if __name__ == '__main__':
   node1 = Ordered_node('C')
   node2 = Ordered_node('C')
   edge = Ordered_edge(node1, node2)
   print(edge.node1_atoms())
