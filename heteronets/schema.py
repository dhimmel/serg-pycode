import networkx


def create_undirected_schema(edges, kind_to_abbrev):
    """edges is an iterable of tuples in the following format:
    (node_kind, edge_kind, node_kind)
    """
    schema = networkx.MultiGraph()
    assert all(len(edge) == 2 for edge in edge_paths)
    node_kinds = set()
    for edge in edges:
        node_kinds.add(edge[0])
        node_kinds.add(edge[2])
    schema.add_nodes_from(node_kinds)
    for edge in edges:
        schema.add_edge(edge[0], edge[2], key=edge[1])
    
    check_abbreviations(schema, kind_to_abbrev)
    schema.graph['paths'] = MetaPaths(self, kind_to_abbrev)
    
    return schema



def check_abbreviations(schema, kind_to_abbrev):
    """Check that the kind to abbreviation dictionary is valid."""
    
    # Check that all kinds have an abbreviation
    kinds_with_abbrev = set(kind_to_abbrev.items())
    kinds = set(schema.nodes() + schema.edges())
    assert kinds <= kinds_with_abbrev
    
    # Check that abbreviations are unique strings of length 1
    abbrevs = kind_to_abbrev.values()
    assert all(isinstance(abbrev, str) for abbrev in abbrevs)
    assert all(len(abbrev) == 1 for abbrev in abbrevs)
    assert len(abbrevs) == len(set(abbrevs))

def extract_metapaths(schema, source, target, max_length,
              exclude_all_source_target_edges = False,
              exclude_edges=set()):
    """
    Calculate paths between a source and target node in the schema.
    exclude_edges is a set of tuples where each tuple is an edge formatted like
    (u, edge, v). Paths that with an edge in exclude_edges are not returned.
    If exclude_all_source_target_edges is True no paths are returned that
    contain an edge between the source and target nodes.
    """
    
    assert isinstance(exclude_edges, set)
    assert all(isinstance(exclude_edge, tuple) for exclude_edge in exclude_edges)
    
    potential_paths = [[source]]
    paths = list()
    for depth in range(max_length):
        current_depth_paths = list()
        for potential_path in potential_paths:
            potential_path_end = potential_path[-1]
            for node, neighbor, key in schema.edges(potential_path_end, keys=True):
                
                # Exclude paths between source and target nodes
                if exclude_all_source_target_edges:
                    if {node, neighbor} == {source, target}:
                        continue

                # Exclude specified paths
                if ((node, key, neighbor) in exclude_edges or
                    (neighbor, key, node) in exclude_edges):
                    continue
                
                path = potential_path + [key, neighbor]
                current_depth_paths.append(path)
                if neighbor == target:
                    path = tuple(path)
                    paths.append(path)
        potential_paths = current_depth_paths
    return paths


class MetaPaths(object):
    
    def __init__(self, schema, kind_to_abbrev):
        self.schema = schema
        
        self.abbrev_to_path = dict()
        #self.paths = set()
        #self.tuple_to_path = dict()
        
        self.kind_to_abbrev = kind_to_abbrev
        self.abbrev_to_kind = {value: key for key, value in self.kind_to_abbrev.iteritems()}
        
    def metapath_from_abbrev(self, abbrev):
        """Return the metapath represented by the supplied abbreviation."""
        metapath = self.abbrev_to_path.get(abbrev)
        if metapath is None:
            kind_tuple = (abbrev_to_kind[x] for x in abbrev)
            metapath = MetaPath(kind_tuple)
        return metapath
    
    def metapath_from_tuple(self, kind_tuple):
        """Return the metapath represented by the supplied tuple of kinds."""
        abbrev = ''.join(self.kind_to_abbrev[kind] for kind in kind_tuple)
        return self.metapath_from_abbrev(abbrev)
        
    def append_to_metapath(self, metapath, edge_kind, node_kind,
                           abbreviated=False, prepend=False):
        """
        Return the metapath representing the supplied metapath appended with
        the supplied edge_kind and node_kind. abbreviated indicates whether the
        supplied node and edge kinds are in their abbreviated form. prepend
        specifies that the node and edge kinds are added to the beginning of
        the metapath.
        """
        if not abbreviated:
            node_kind = self.abbrev_to_kind[node_kind]
            edge_kind = self.abbrev_to_kind[edge_kind]
        if prepend:
            abbrev = node_kind + edge_kind + metapath.abbrev
        else:
            abbrev = metapath.abbrev + edge_kind + node_kind
        return self.metapath_from_abbrev(abbrev)
        

class MetaPath(object):
    
    def __init__(self, metapaths, kind_tuple):
        self.tuple_ = kind_tuple
        self.abbrev = ''.join(metapaths.kind_to_abbrev[kind] for kind in kind_tuple)
        assert self.abbrev not in metapaths.abbrev_to_path
        metapaths.abbrev_to_path[self.abbrev] = self
        
    def __len__(self):
        return len(self.abbrev) / 2
    
    def __hash__(self):
        return hash(self.abbrev)
    
    def __eq__(self, other):
        return self.abbrev == other.abbrev
    
    def __str__(self):
        return self.abbrev

if __name__ == '__main__':
    
    """
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('gene', 'gene', 'function')]
    edge_tuples = [('drug', 'gene', 'target'),
                   ('gene', 'disease', 'risk'),
                   ('drug', 'disease', 'indication')]

    edge_tuples = [('drug', 'gene', 'upregulates'),
                   ('drug', 'gene', 'downregulates'),
                   ('disease', 'gene', 'upregulates'),
                   ('disease', 'gene', 'downregulates'),
                   ('gene', 'gene', 'function')]
    """
    