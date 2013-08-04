import logging
import networkx


def create_undirected_schema(edges, kind_to_abbrev=None):
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
    
    if not kind_to_abbrev:
        kind_to_abbrev = create_abbreviations(node_kinds, edges)
    else:
        check_abbreviations(schema, kind_to_abbrev)
    MetaPath.kind_to_abbrev = kind_to_abbrev    
    log_undirected_schema(schema)
    return schema

def log_undirected_schema(schema):
    """ """
    kind_to_abbrev = MetaPath.kind_to_abbrev
    
    lines = []
    lines.append('----------Undirected Graph Schema------------')
    lines.append('Nodes (abbreviation - kind):')
    for node in schema.nodes():
        abbrev = kind_to_abbrev[node]
        lines.append(abbrev + ' - ' + node)
    
    lines.append('Edges (abbreviation - kind_tuple):')
    for edge in schema.edges(keys=True):
        abbrev = kind_to_abbrev[edge[2]]
        lines.append(abbrev + ' - ' + str(edge))
    
    msg = '\n'.join(lines)
    logging.info(msg)
    

def get_duplicates(iterable):
    """Return a set of the elements which appear multiple times in iterable."""
    seen, duplicates = set(), set()
    for elem in iterable:
        if elem in seen:
            duplicates.add(elem)
        else:
            seen.add(elem)
    return duplicates

def find_abbrevs(kinds):
    """For a list of strings (kinds), find the shortest unique abbreviation."""
    kind_to_abbrev = {kind: kind[0] for kind in kinds}
    duplicates = get_duplicates(kind_to_abbrev.values())
    while duplicates:
        for kind, abbrev in kind_to_abbrev.items():
            if abbrev in duplicates:
                abbrev += kind[len(abbrev)]
                kind_to_abbrev[kind] = abbrev
        duplicates = get_duplicates(kind_to_abbrev.values())
    return kind_to_abbrev
    
def create_abbreviations(node_kinds, edge_kinds):
    """Creates abbreviations for node and edge kinds."""
    kind_to_abbrev = find_abbrevs(node_kinds)
    kind_to_abbrev = {kind: abbrev.upper()
                      for kind, abbrev in kind_to_abbrev.items()}
    
    edge_set_to_keys = dict()
    for edge in edge_kinds:
        key = frozenset(map(str.lower, edge[:2]))
        value = edge[2]
        edge_set_to_keys.setdefault(key, list()).append(value)
    
    for edge_set, keys in edge_set_to_keys.items():
        key_to_abbrev = find_abbrevs(keys)
        for key, abbrev in key_to_abbrev.items():
            previous_abbrev = kind_to_abbrev.get(key)
            if previous_abbrev and len(abbrev) <= len(previous_abbrev):
                continue
            kind_to_abbrev[key] = abbrev
    
    return kind_to_abbrev

    

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
    kind_to_abbrev = MetaPath.kind_to_abbrev
    for kind, abbrev in kind_to_abbrev.items():
        print kind + ': ' + abbrev

def check_abbreviations(schema, kind_to_abbrev):
    """Check that the kind to abbreviation dictionary is valid."""
    
    # Check that all kinds have an abbreviation
    kinds_with_abbrev = set(kind_to_abbrev)
    kinds = set(schema.nodes()) | {key for n1, n2, key in schema.edges(keys=True)}
    assert kinds <= kinds_with_abbrev
    
    # Check that abbreviations are strings
    abbrevs = kind_to_abbrev.values()
    assert all(isinstance(abbrev, str) for abbrev in abbrevs)

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
            
    potential_paths = [MetaPath((source, ))]
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
                
                path = potential_path.append(key, neighbor)
                
                current_depth_paths.append(path)
                if neighbor == target:
                    paths.append(path)
        potential_paths = current_depth_paths
    return paths


class MetaPath(object):
    
    tuple_to_metapath = dict()
    kind_to_abbrev = dict()

    def __new__(cls, kind_tuple):
        """
        If a MetaPath object representing kind_tuple exists, return that
        object. Otherwise create and return a new MetaPath object.
        """
        assert isinstance(kind_tuple, tuple)
        
        if kind_tuple in MetaPath.tuple_to_metapath:
            return MetaPath.tuple_to_metapath[kind_tuple]

        new = object.__new__(cls)
        new.tuple_ = kind_tuple
        new.hash_ = hash(kind_tuple)
        #new.abbrev = new.get_abbrev()
        MetaPath.tuple_to_metapath[kind_tuple] = new
        return new

    def __getnewargs__(self):
        """Needed for pickling."""
        return (self.tuple_, )
        
    @staticmethod
    def edges_in_metapaths(metapaths):
        """Returns a set of edges appearing in any of the metapaths"""
        edges = set()
        for metapath in metapaths:
            edges |= set(metapath.get_edges())
        return edges    

    @staticmethod
    def shortcuts_for_metapaths(metapaths, shortcut_length):
        """Compute desired shortcuts for faster computation of metapaths."""
        shortcuts = set()
        for metapath in metapaths:
            shortcuts |= metapath.shortcuts(shortcut_length)
        return shortcuts

    def shortcuts(self, length):
        """Get shortcuts for self"""
        shorcuts = set()
        num_nodes = len(self)
        depth = 0
        while depth + length <= num_nodes:
            start = depth * 2
            end = start + length * 2 + 1
            shortcut_tuple = self[start: end]
            shortcut = MetaPath(shortcut_tuple)
            shorcuts.add(shortcut)
            depth += length
        return shorcuts
    
    def get_edges(self):
        """Return a the list of edges composing the metapath."""
        if hasattr(self, 'edges'):
            return self.edges
        
        edges = list()
        for depth in xrange(len(self)):
            edge = self[depth * 2 : depth * 2 + 3]
            edge = edge[0], edge[2], edge[1]
            edges.append(edge)
        
        self.edges = edges
        return self.edges
    
    def get_abbrev(self):
        """Return the abbreviation representing the metapath"""
        if not hasattr(self, 'abbrev'):        
            abbrev_generator = (MetaPath.kind_to_abbrev[kind] for kind in self)
            self.abbrev = str.join('', abbrev_generator)
        return self.abbrev
    
    def startswith(self, other):
        return self[: len(other.tuple_)] == other.tuple_

    def endswith(self, other):
        return self[-len(other.tuple_): ] == other.tuple_
    
    def start(self):
        """Return the node_kind for the start node of the metapath."""
        return self[0]

    def end(self):
        """Same as start except returns the node_kind for the end node."""
        return self[-1]
    
    def append(self, edge_kind, node_kind, prepend=False):
        """
        Return the metapath representing the supplied metapath appended with
        the supplied edge_kind and node_kind. prepend specifies that the node
        and edge kinds are added to the beginning of the metapath.
        """
        if prepend:
            kind_tuple = (node_kind, edge_kind) + self.tuple_
        else:
            kind_tuple = self.tuple_ + (edge_kind, node_kind)
        return MetaPath(kind_tuple)
    
    def split_by_index(self, index, reverse_head = False):
        """Split the metapath into two metapaths. Index is included as the 
        end node kind of head and the start node kind of tail.
        """
        head = MetaPath(self[: index * 2 + 1])
        tail = MetaPath(self[index * 2: ])
        if reverse_head:
            head = head.reverse()
        return head, tail
    
    def split_by_kind(self, kind, reverse_head = False):
        """Return a list of (head, tail) tuples where the tuples are self split
        at each occurance of node kind.
        """
        nodes = self.tuple_[::2]
        split_tuples = [self.split_by_index(i, reverse_head)
                for i in range(len(nodes)) if nodes[i] == kind]
        return split_tuples
    
    def reverse(self):
        return MetaPath(tuple(reversed(self.tuple_)))
    
    def __reversed__(self):
        """For creating a reversed MetaPath object use MetaPath.reverse().
        Iterators don't provide a meaningful representation.
        """
        raise Exception
    
    def __getitem__(self, index):
        return self.tuple_[index]
    
    def __iter__(self):
        return iter(self.tuple_)
    
    def __len__(self):
        return len(self.tuple_) / 2
    
    def __hash__(self):
        return self.hash_
    
    def __eq__(self, other):
        return self.tuple_ == other.tuple_
    
    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return self.get_abbrev()

if __name__ == '__main__':
    
    
    edges = [('drug', 'gene', 'target'),
             ('gene', 'disease', 'risk'),
             ('drug', 'disease', 'indication'),
             ('drug', 'disease', 'increases'),
             ('gene', 'gene', 'function')]
    
    kind_to_abbrev = {'drug': 'C', 'disease': 'D', 'gene': 'G',
                      'risk': 'r', 'indication': 'i', 'target': 't', 'function': 'f', 'increases': '^'}
    
    edges = [('compound', 'gene', 'up-regulation'),
             ('compound', 'gene', 'down-regulation'),
             ('disease', 'gene', 'up-regulation'),
             ('disease', 'gene', 'down-regulation'),
             ('disease', 'gene', 'regulation'),
             ('compound', 'gene', 'regulation'),
             ('gene', 'gene', 'interaction'),
             ('gene', 'tissue', 'expression'),
             ('compound', 'side effect', 'causation'),
             ('disease', 'disease', 'comorbidity'),
             ('compound', 'compound', 'similarity'),
             ('disease', 'gene', 'association'),
             ('disease', 'compound', 'indication'),
             ]
    
    schema = create_undirected_schema(edges)
    
    #print_schema(schema)
    
    source = 'compound'
    target = 'disease'
    max_length = 3
    
    print 'Metapaths'
    metapaths = extract_metapaths(schema, source, target, 3,
        exclude_edges={('disease', 'gene', 'regulation'), ('compound', 'gene', 'regulation')})
    #print metapaths    
    metapaths = extract_metapaths(schema, source, target, 4,
        exclude_edges={('disease', 'gene', 'up-regulation'), ('disease', 'gene', 'down-regulation'),
                       ('compound', 'gene', 'up-regulation'), ('compound', 'gene', 'down-regulation')})
    print metapaths
    print len(metapaths)    

    #print 'Shortcuts'
    #shortcuts = MetaPath.shortcuts_for_metapaths(metapaths, 2)
    #print shortcuts
    
    #print 'Metapaths Excluding source to target edges'
    #metapaths = extract_metapaths(schema, source, target, max_length, exclude_all_source_target_edges=True)
    #print metapaths

    #print 'Metapaths Excluding (gene, function, gene) edges'
    #metapaths = extract_metapaths(schema, source, target, max_length, exclude_edges={('gene', 'gene', 'function')})
    #print metapaths

    #print 'Metapaths Excluding (disease, indication, drug) edges'
    #metapaths = extract_metapaths(schema, source, target, max_length, exclude_edges={('disease', 'drug', 'indication')})
    
    
    # ASHG Network
    edges = [('disease', 'gene', 'up-regulation'),
             ('disease', 'gene', 'down-regulation'),
             ('disease', 'gene', 'association'),
             ('disease', 'disease', 'comorbidity'),
             ('gene', 'gene', 'interaction')
             ]

    edges = [('disease', 'gene', 'regulation'),
             ('disease', 'gene', 'association'),
             ('disease', 'disease', 'comorbidity'),
             ('gene', 'gene', 'interaction')
             ]
    
    schema = create_undirected_schema(edges)
    metapaths = extract_metapaths(schema, 'disease', 'gene', 3)
    #print_schema(schema)
    #print metapaths
    