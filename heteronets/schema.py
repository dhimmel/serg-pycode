import networkx


def create_undirected_schema(edges, kind_to_abbrev):
    """
    Returns a networkx.MultiGraph graph representing the schema describing a
    heterogeneous information network. 
    
    edges is an iterable of tuples in the following format:
    (node_kind, node_kind, edge_kind)
    
    kind_to_abbrev is a dictionary where the keys refer to kinds input in edges
    and the values refer to the unique one character string abbreviation for
    that edge.
    """
    schema = networkx.MultiGraph()
    assert all(len(edge) == 3 for edge in edges)
    node_kinds = set()
    for edge in edges:
        node_kinds.add(edge[0])
        node_kinds.add(edge[1])
    schema.add_nodes_from(node_kinds)
    for edge in edges:
        schema.add_edge(edge[0], edge[1], key=edge[2])
    
    check_abbreviations(schema, kind_to_abbrev)
    schema.graph['paths'] = MetaPaths(schema, kind_to_abbrev)
    
    return schema

def print_schema(schema):
    """Print schema information to console."""
    
    print 'Schema Nodes:'
    for node in schema.nodes():
        print node

    print 'Schema Edges'
    edges = schema.edges(keys=True)
    edges.sort()
    for edge in edges:
        print edge
    
    print 'Kind Abbreviations'
    kind_to_abbrev = schema.graph['paths'].kind_to_abbrev
    for kind, abbrev in kind_to_abbrev.items():
        print kind + ': ' + abbrev

def check_abbreviations(schema, kind_to_abbrev):
    """Check that the kind to abbreviation dictionary is valid."""
    
    # Check that all kinds have an abbreviation
    kinds_with_abbrev = set(kind_to_abbrev)
    kinds = set(schema.nodes()) | {key for n1, n2, key in schema.edges(keys=True)}
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
    (u, v, edge). Paths that with an edge in exclude_edges are not returned.
    If exclude_all_source_target_edges is True no paths are returned that
    contain an edge between the source and target nodes.
    """
    
    assert isinstance(exclude_edges, set)
    assert all(isinstance(exclude_edge, tuple) for exclude_edge in exclude_edges)
    
    paths_obj = schema.graph['paths']
        
    potential_paths = [paths_obj.metapath_from_tuple((source, ))]
    paths = list()
    for depth in range(max_length):
        current_depth_paths = list()
        for potential_path in potential_paths:
            for node, neighbor, key in schema.edges(potential_path.end(), keys=True):
                
                # Exclude paths between source and target nodes
                if (exclude_all_source_target_edges and
                    {node, neighbor} == {source, target}):
                    continue
                # Exclude specified paths
                if ((node, neighbor, key) in exclude_edges or
                    (neighbor, node, key) in exclude_edges):
                    continue
                
                path = paths_obj.append_to_metapath(potential_path, key, neighbor)
                
                current_depth_paths.append(path)
                if neighbor == target:
                    paths.append(path)
        potential_paths = current_depth_paths
    return paths

def shortcuts_for_metapaths(schema, metapaths, shortcut_length):
    """Compute desired shortcuts for faster computation of metapaths."""
    shorcuts = set()
    for metapath in metapaths:
        num_nodes = len(metapath)
        depth = 0
        while depth + shortcut_length <= num_nodes:
            start = depth * 2
            end = start + shortcut_length * 2 + 1
            shortcut = metapath.tuple_[start:end]
            shortcut = schema.graph['paths'].metapath_from_tuple(shortcut)
            shorcuts.add(shortcut)
            depth += shortcut_length
    return shorcuts

def edges_in_metapaths(metapaths):
    """Returns a set of edges appearing in any of the metapaths"""
    edges = set()
    for metapath in metapaths:
        edges |= metapath.get_edges()
    return edges

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
            kind_tuple = tuple(self.abbrev_to_kind[x] for x in abbrev)
            metapath = MetaPath(self, kind_tuple)
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
            node_kind = self.kind_to_abbrev[node_kind]
            edge_kind = self.kind_to_abbrev[edge_kind]
        if prepend:
            abbrev = node_kind + edge_kind + metapath.abbrev
        else:
            abbrev = metapath.abbrev + edge_kind + node_kind
        return self.metapath_from_abbrev(abbrev)
        

class MetaPath(object):
    
    def __init__(self, metapaths, kind_tuple):
        self.metapaths = metapaths
        self.tuple_ = kind_tuple
        self.abbrev = ''.join(metapaths.kind_to_abbrev[kind] for kind in kind_tuple)
        assert self.abbrev not in metapaths.abbrev_to_path
        metapaths.abbrev_to_path[self.abbrev] = self
    
    def start(self):
        """Return the node_kind for the start node of the metapath."""
        return self.tuple_[0]

    def end(self):
        """Same as start except returns the node_kind for the end node."""
        return self.tuple_[-1]
    
    def get_edges(self):
        """Return a the set of edges composing the metapath."""
        if hasattr(self, 'edges'):
            return self.edges
        
        edges = set()
        for depth in xrange(len(self)):
            edge = self.tuple_[depth * 2 : depth * 2 + 3]
            edge = edge[0], edge[2], edge[1]
            edges.add(edge)
        
        self.edges = edges
        return self.edges
    
    def split_by_index(self, index, reverse_head = False):
        """Split the metapath into two metapaths. Index is included as the 
        end node kind of head and the start node kind of tail.
        """
        head = self.metapaths.metapath_from_tuple(self.tuple_[: index * 2 + 1])
        tail = self.metapaths.metapath_from_tuple(self.tuple_[index * 2: ])
        if reverse_head:
            head = reversed(head)
        return head, tail
    
    def split_by_kind(self, kind, reverse_head = False):
        """Return a list of (head, tail) tuples where the tuples are self split
        at each occurance of node kind.
        """
        nodes = self.tuple_[::2]
        split_tuples = [self.split_by_index(i, reverse_head)
                for i in range(len(nodes)) if nodes[i] == kind]
        return split_tuples
    
    def __reversed__(self):
        return self.metapaths.metapath_from_tuple(reversed(self.tuple_))
    
    def __len__(self):
        return len(self.tuple_) / 2
    
    def __hash__(self):
        return hash(self.abbrev)
    
    def __eq__(self, other):
        return self.abbrev == other.abbrev
    
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.abbrev

if __name__ == '__main__':
    
    edges = [('drug', 'gene', 'target'),
             ('gene', 'disease', 'risk'),
             ('drug', 'disease', 'indication'),
             ('drug', 'disease', 'increases'),
             ('gene', 'gene', 'function')]
    
    kind_to_abbrev = {'drug': 'C', 'disease': 'D', 'gene': 'G',
                      'risk': 'r', 'indication': 'i', 'target': 't', 'function': 'f', 'increases': '^'}
    
    
    schema = create_undirected_schema(edges, kind_to_abbrev)
    
    print_schema(schema)
    
    source = 'drug'
    target = 'disease'
    max_length = 5
    
    print 'Metapaths'
    metapaths = extract_metapaths(schema, source, target, max_length)
    shortcuts = shortcuts_for_metapaths(schema, metapaths, 2)
    print metapaths
    print 'Shortcuts'
    print shortcuts
    
    print 'Metapaths Excluding source to target edges'
    metapaths = extract_metapaths(schema, source, target, max_length, exclude_all_source_target_edges=True)
    print metapaths

    print 'Metapaths Excluding (gene, function, gene) edges'
    metapaths = extract_metapaths(schema, source, target, max_length, exclude_edges={('gene', 'gene', 'function')})
    print metapaths

    print 'Metapaths Excluding (disease, indication, drug) edges'
    metapaths = extract_metapaths(schema, source, target, max_length, exclude_edges={('disease', 'drug', 'indication')})
    print metapaths
    
    
    print metapaths[2].get_edges()
    print metapaths[2].split(0)
    print metapaths[2].split(1)
    
    """
    edge_tuples = [('drug', 'gene', 'upregulates'),
                   ('drug', 'gene', 'downregulates'),
                   ('disease', 'gene', 'upregulates'),
                   ('disease', 'gene', 'downregulates'),
                   ('gene', 'gene', 'function')]
    """
    
    
    