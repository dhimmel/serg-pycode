import networkx
import matplotlib.pyplot

class Schema(object):
    
    def __init__(self, g, node_kinds, edge_tuples):
        """Create a graph schema representing node and edge kinds. Meant for
        either directed or undirected multigraphs."""
        self.g = g
        g.graph['kind'] = 'schema'
        g.add_nodes_from(node_kinds)
        self.node_kinds = node_kinds
        self.edge_kinds = set()
        for source, target, kind in edge_tuples:
            self.edge_kinds.add(kind)
            g.add_edge(source, target, key=kind)
        
    def plot(self):
        networkx.draw_circular(self.g)
        matplotlib.pyplot.show()
    
    
class UndirectedSchema(Schema):
    
    def __init__(self, node_kinds, edge_tuples):
        """Create a graph schema representing node and edge kinds for an
        undirected graph."""
        g = networkx.MultiGraph()
        super(UndirectedSchema, self).__init__(g, node_kinds, edge_tuples)
        
    def metapaths(self, source, target, cutoff=3):
        paths = networkx.all_simple_paths(self.g, source, target, cutoff)
        potential_paths = [[source]]
        paths = list()
        for depth in range(cutoff):
            current_depth_paths = list()
            for potential_path in potential_paths:
                for node, neighbor, key in self.g.edges(potential_path[-1], keys=True):
                    path = potential_path + [key, neighbor]
                    current_depth_paths.append(path)
                    if neighbor == target:
                        path = tuple(path)
                        paths.append(path)
            potential_paths = current_depth_paths        
        return paths



if __name__ == '__main__':
    node_kinds = {'drug', 'disease', 'gene'}
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('gene', 'gene', 'function')]
    omicnet_schema = UndirectedSchema(node_kinds, edge_tuples)
    #print omicnet_schema.g.edges('gene')
    print list(omicnet_schema.metapaths('drug', 'disease'))