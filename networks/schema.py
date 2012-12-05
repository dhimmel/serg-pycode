import networkx
import matplotlib.pyplot

class Schema(object):
    
    def __init__(self, g, node_kinds, edge_tuples):
        """Create a graph schema representing node and edge kinds. Meant for
        either directed or undirected multigraphs."""
        self.edge_tuples = edge_tuples
        self.g = g
        g.graph['kind'] = 'schema'
        g.add_nodes_from(node_kinds)
        self.node_kinds = node_kinds
        self.edge_kinds = set()
        for source, target, kind in edge_tuples:
            self.edge_kinds.add(kind)
            assert source in g and target in g
            g.add_edge(source, target, key=kind)
        
    def plot(self):
        networkx.draw_circular(self.g)
        matplotlib.pyplot.show()
    
    def __str__(self):
        edge_list = list()
        for edge_tuple in self.edge_tuples:
            edge = [edge_tuple[0], edge_tuple[2], edge_tuple[1]]
            edge_list.append(self.path_as_str(edge))
        return 'Nodes: ' + ', '.join(self.node_kinds) + '\nEdges:\n' + '\n'.join(edge_list)
    
    @staticmethod
    def path_as_str(path):
        path = list(path)
        for i in range(len(path)):
            if i % 2:
                path[i] = ' <--' +  path[i] + '--> '
        return ''.join(path)
    
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


def path_count(g, source, target, metapath):
    """ """
    assert g.node[source]['kind'] == metapath[0]
    assert g.node[target]['kind'] == metapath[-1]
    metapath.pop(0)
    paths = [[source]]
    while metapath:
        edge_key = metapath.pop(0)
        node_kind = metapath.pop(0)
        current_depth_paths = list()
        while paths:
            preceeding_path = current_depth_paths.pop()
            node = preceeding_path[-1]
            
            for node, neighbor, key in g.edges(node, keys=True):
                neighbor_kind = g.node[neighbor]['kind']
                if key == edge_key and neighbor_kind == node_kind:
                    path = preceeding_path + [edge_key, neighbor]
                    current_depth_paths.append(path)
        paths = current_depth_paths
    paths_from_source = len(paths)
    paths = filter(lambda path: path[-1] == target, paths)
    return len(paths), paths

if __name__ == '__main__':
    node_kinds = {'drug', 'disease', 'gene'}
    """
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('gene', 'gene', 'function')]
    """
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('drug', 'disease', 'indication')]

    
    """
    edge_tuples = [('drug', 'gene', 'upregulates'),
                   ('drug', 'gene', 'downregulates'),
                   ('disease', 'gene', 'upregulates'),
                   ('disease', 'gene', 'downregulates'),
                   ('gene', 'gene', 'function')]
    """
    
    omicnet_schema = UndirectedSchema(node_kinds, edge_tuples)
    #print omicnet_schema.g.edges('gene')
    metapaths = omicnet_schema.metapaths('drug', 'disease', 3)
    metapath_strs = [Schema.path_as_str(path) for path in metapaths]
    print omicnet_schema
    print 'Metapaths with length <= 3:'
    for metapath in metapaths:
        print Schema.path_as_str(metapath)
        
        
