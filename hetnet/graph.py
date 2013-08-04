# daniel.himmelstein@gmail.com

class Schema(object):
    
    inverse_direction = {'forward': 'backward',
                         'backward': 'forward',
                         'both': 'both'}
    
    direction_to_abbrev = {'forward': '>', 'backward': '<', 'both': '-'}
    
    def __init__(self, edge_kinds, kind_to_abbrev=None):
        """edges is a list of tuples in the format: (source, target, kind, direction)"""
        self.edge_kinds = set()
        self.node_kinds = set()

        for edge_kind in edge_kinds:
            source_kind, target_kind, kind, direction = edge_kind

            for node_kind in [source_kind, target_kind]:
                self.node_kinds.add(node_kind)

            inverse_edge_kind = target_kind, source_kind, kind, Schema.inverse_direction[direction]
            for ek in [edge_kind, inverse_edge_kind]:
                self.edge_kinds.add(ek)

        self.edges_to_metapath = dict()
        
        if kind_to_abbrev is None:
            kind_to_abbrev = self.create_abbreviations()
        self.kind_to_abbrev = kind_to_abbrev
    
    def get_metapath(self, edges):
        """Retreive or construct metapath based on edges"""
        try:
            return self.edges_to_metapath[edges]
        except KeyError:
            metapath = MetaPath(self, edges)
            self.edges_to_metapath[edges] = metapath
            return metapath
    
    def extract_metapaths(self, source, target, max_length):
        metapaths = list()
        metapaths_of_depth = list()
        for depth in range(max_length):
            if depth == 0:
                metapaths_of_depth = [MetaPath(self, (e, )) for e in self.edge_kinds if e[0] == source]
            else:
                temp_paths = list()
                for metapath in metapaths_of_depth:
                    node_kind = metapath.get_nodes()[-1]
                    temp_paths.extend(metapath.append(e) for e in self.edge_kinds if e[0] == node_kind)
                metapaths_of_depth = temp_paths
            
            for metapath in metapaths_of_depth:
                if target == metapath.get_nodes()[-1]:
                    metapaths.append(metapath)
        
        return metapaths

    @staticmethod
    def get_duplicates(iterable):
        """Return a set of the elements which appear multiple times in iterable."""
        seen, duplicates = set(), set()
        for elem in iterable:
            if elem in seen:
                duplicates.add(elem)
            else:
                seen.add(elem)
        return duplicates
    
    @staticmethod
    def find_abbrevs(kinds):
        """For a list of strings (kinds), find the shortest unique abbreviation."""
        kind_to_abbrev = {kind: kind[0] for kind in kinds}
        duplicates = Schema.get_duplicates(kind_to_abbrev.values())
        while duplicates:
            for kind, abbrev in kind_to_abbrev.items():
                if abbrev in duplicates:
                    abbrev += kind[len(abbrev)]
                    kind_to_abbrev[kind] = abbrev
            duplicates = Schema.get_duplicates(kind_to_abbrev.values())
        return kind_to_abbrev

    def create_abbreviations(self):
        """Creates abbreviations for node and edge kinds."""
        kind_to_abbrev = Schema.find_abbrevs(self.node_kinds)
        kind_to_abbrev = {kind: abbrev.upper()
                          for kind, abbrev in kind_to_abbrev.items()}
        
        edge_set_to_keys = dict()
        for edge in self.edge_kinds:
            key = frozenset(map(str.lower, edge[:2]))
            value = edge[2]
            edge_set_to_keys.setdefault(key, list()).append(value)
        
        for edge_set, keys in edge_set_to_keys.items():
            key_to_abbrev = Schema.find_abbrevs(keys)
            for key, abbrev in key_to_abbrev.items():
                previous_abbrev = kind_to_abbrev.get(key)
                if previous_abbrev and len(abbrev) <= len(previous_abbrev):
                    continue
                kind_to_abbrev[key] = abbrev
        
        return kind_to_abbrev

    def get_node_dicts(self, node_kind):
        edge_kind_filter = lambda d: (ek[2] for ek in self.edge_kinds
                                      if ek[0] == node_kind and ek[3] == d)
        incoming = {key: set() for key in edge_kind_filter('backward')}
        outgoing = {key: set() for key in edge_kind_filter('forward')}
        incident = {key: set() for key in edge_kind_filter('both')}
        return incoming, outgoing, incident

class Path(object):
    
    def __init__(self, edges):
        self.edges = edges
    
    def __hash__(self):
        return hash(self.edges)

    def __eq__(self, other):
        return self.edges == other.edges

    def get_nodes(self):
        nodes = [edge[0] for edge in self.edges]
        nodes.append(self.edges[-1][1])
        return tuple(nodes)

    def get_edge_kinds(self):
        return tuple(edge[2] for edge in self.edges)

    def get_directions(self):
        return tuple(edge[3] for edge in self.edges)

    def __repr__(self):
        return str(self.edges)

        
class MetaPath(Path):
    
    def __init__(self, schema, edge_kinds):
        """schema is the Schema object defining the metapath. kind is a tuple
        of the node and edge kinds composing the metapath."""
        super(MetaPath, self).__init__(edge_kinds)
        self.edge_kinds = self.edges
        self.schema = schema
        for edge_kind in edge_kinds:
            assert edge_kind in schema.edge_kinds

    def append(self, edge_kind):
        """Returns a new MetaPath with the edge appended."""
        edge_kinds = self.edge_kinds + (edge_kind, )
        return self.schema.get_metapath(edge_kinds)
    
    def __repr__(self):
        if hasattr(self, 'abbreviation'):
            return self.abbreviation
        kind_to_abbrev = self.schema.kind_to_abbrev
        nodes = [kind_to_abbrev[x] for x in self.get_nodes()]
        edges = [kind_to_abbrev[x] for x in self.get_edge_kinds()]
        directions = [Schema.direction_to_abbrev[x] for x in self.get_directions()]
        abbreviation = ''
        for i in range(len(edges)):
            direction = directions[i]
            abbreviation += nodes[i] + direction + edges[i] + direction
        abbreviation += nodes[-1]
        self.abbreviation = abbreviation
        return abbreviation
            
    
class Graph(object):
    
    def __init__(self, schema, data=dict()):
        """Graph class"""
        self.data = data
        self.node_dict = dict()
        self.edge_dict = dict()
        self.schema = schema

    def get_node(self, node_id):
        """ """
        return self.node_dict[node_id]
        
    def get_edge(self, source, target, kind, direction):
        """ """
        edge_id = source, target, kind, direction
        return self.edge_dict[edge_id]
        
    def get_nodes(self):
        """iterator over nodes"""
        return self.node_dict.itervalues()
    
    def get_edges(self, include_inverts=False):
        """iterator over edges"""
        if include_inverts:
            return self.edge_dict.itervalues()
        else:
            return (edge for edge in self.edge_dict.itervalues() if edge.inverted == False)
    
    def get_paths(self, source, target, metapath):
        """ """
        paths = [Path(tuple())]
        for i in range(len(metapath)):
            source_kind, target_kind, edge_kind, direction = metapath[i]
            for path in paths:
                nodes = path.get_nodes()
                node = nodes[-1]
                attribute = Node.direction_to_attr[direction]
                nodes[-1].getattr
    
    def add_node(self, node_id, kind, data=dict()):
        """ """
        node = Node(self, node_id, kind, data)
        return node
    
    def add_edge(self, source_id, target_id, kind, direction, data=dict()):
        """ """
        source = self.node_dict[source_id]
        target = self.node_dict[target_id]
        edge = Edge(self, source, target, kind, direction, data)
        inverse = edge.invert()
        return edge, inverse

    def remove_unconnected_nodes(self):
        unconnected = [node for node in self.get_nodes() if node.degree() == 0]
        for node in unconnected:
            node.remove()

class Node(object):
    
    direction_to_attr = {'both': 'incident', 'forward': 'outgoing', 'backward': 'incoming'}

    def __init__(self, graph, node_id, kind, data):
        """Graph class"""
        self.__dict__.update(locals())
        assert kind in graph.schema.node_kinds
        assert node_id not in graph.node_dict
        graph.node_dict[node_id] = self

        node_dicts = graph.schema.get_node_dicts(kind)
        self.incoming, self.outgoing, self.incident = node_dicts
    
    def all_edges(self):
        all_edges = list()
        for kind_to_edges in (self.incoming, self.outgoing, self.incident):
            for edges in kind_to_edges.values():
                all_edges.extend(edges)
        return all_edges
    
    def degree(self):
        return len(self.all_edges())
    
    def __hash__(self):
        return hash(self.node_id)
    
    def __eq__(self, other):
        return self.node_id == other.node_id
    
    def __repr__(self):
        return self.node_id
    
    def remove(self):
        for edge in self.all_edges():
            edge.remove()
        del self.graph.node_dict[self.node_id]
    
class Edge(object):
    
    
    def __init__(self, graph, source, target, kind, direction, data):
        """Graph class"""
        self.__dict__.update(locals())
        assert self.get_edge_kind() in graph.schema.edge_kinds
        edge_id = self.get_edge_id()
        assert edge_id not in graph.edge_dict
        graph.edge_dict[edge_id] = self
        getattr(source, Node.direction_to_attr[direction])[kind].add(self)

        
    def get_edge_kind(self):
        return self.source.kind, self.target.kind, self.kind, self.direction
    
    def get_edge_id(self):
        return self.source.node_id, self.target.node_id, self.kind, self.direction
    
    def invert(self):
        inverse = Edge(self.graph, self.target, self.source, self.kind,
                    Schema.inverse_direction[self.direction], self.data)
        inverse.inverted = True
        self.inverted = False
        self.inverse = inverse
        inverse.inverse = self
        return inverse
    
    def remove(self):
        for edge in (self, self.inverse):
            getattr(edge.source, Node.direction_to_attr[edge.direction])[edge.kind].remove(edge)
            del edge.inverse
            del edge.graph.edge_dict[edge.get_edge_id()]
    
    def __hash__(self):
        return hash(self.get_edge_id())
    
    def __eq__(self, other):
        return self.get_edge_id() == other.get_edge_id()
    
    def __repr__(self):
        direction_abbrev = Schema.direction_to_abbrev[self.direction]
        return '%s %s %s %s %s' % (self.source, direction_abbrev, self.kind, direction_abbrev, self.target)

if __name__ == '__main__':
    """ """
    schema_edges = [('gene', 'disease', 'association', 'both'),
             ('gene', 'gene', 'function', 'both'),
             ('gene', 'tissue', 'expression', 'both'),
             ('disease', 'tissue', 'pathology', 'both'),
             ('gene', 'gene', 'transcription', 'forward')]

    schema = Schema(schema_edges)
    graph = Graph(schema)
    graph.add_node('IL17', 'gene')
    graph.add_node('BRCA1', 'gene')
    graph.add_node('ANAL1', 'gene')
    graph.add_node('MS', 'disease')
    graph.add_node('brain', 'tissue')
    
    graph.add_edge('MS', 'IL17', 'association', 'both')
    graph.add_edge('MS', 'brain', 'pathology', 'both')
    graph.add_edge('IL17', 'brain', 'expression', 'both')
    graph.add_edge('IL17', 'BRCA1', 'transcription', 'backward')

    for node in graph.get_nodes():
        print node, node.degree()

    

    #metapath_edges = (('gene', 'disease', 'association', 'both'),)
    #schema.metapath(metapath_edges)
    #print schema.extract_metapaths('gene', 'disease', 3)
    
    