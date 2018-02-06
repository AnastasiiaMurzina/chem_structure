class Ordered_node():
    order_edges = []
    element = ''

    def __init__(self, name, order_edges = []):
        self.element = name
        self.order_edges = order_edges

    # def __init__(self, count_max_bounds: object, order_edges: object) -> object:
    #     self.count_max_bounds = count_max_bounds
    #     self.order_edges = order_edges


if __name__ == '__main__':
    node = Ordered_node('C')
    node1 = Ordered_node('H')
    node2 = Ordered_node('H')
    node.order_edges.append(node1)
    node.order_edges.append(node2)
    from Ordered_edge import Ordered_edge
    edge = Ordered_edge(node, node1)
    # for i in edge.node1.order_edges:
    #     print(i.element)
    # print(edge.node1.order_edges)
