class BasePath(object):
    
    def __init__(self, edges):
        assert isinstance(edges, tuple)
        self.edges = edges
    
    def source(self):
        return self[0].source

    def target(self):
        return self[-1].target
    
    def get_nodes(self):
        nodes = tuple(edge.source for edge in self)
        nodes = nodes + (self.target(), )
        return nodes
    
    def __iter__(self):
        return iter(self.edges)

    def __getitem__(self, key):
        return self.edges[key]

    def __len__(self):
        return len(self.edges)

    def __add__(self, other):
        assert type(self) is type(other) 
        return type(self)(self.edges + other.edges)

class MetaGraph(object):
    
    def __init__(self):
        """ """
        self.node_dict = dict()
        self.edge_dict = dict()
        self.path_dict = dict()
    
    
    @staticmethod
    def from_edge_tuples(metaedge_tuples):
        metagraph = MetaGraph()
        node_kinds = set()
        for source_kind, target_kind, kind, direction in metaedge_tuples:
            node_kinds.add(source_kind)
            node_kinds.add(target_kind)
        for kind in node_kinds:
            metagraph.add_node(kind)
        for edge_tuple in metaedge_tuples:
            metagraph.add_edge(edge_tuple)
        
        metagraph.create_abbreviations()
        return metagraph

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
        duplicates = MetaGraph.get_duplicates(kind_to_abbrev.values())
        while duplicates:
            for kind, abbrev in kind_to_abbrev.items():
                if abbrev in duplicates:
                    abbrev += kind[len(abbrev)]
                    kind_to_abbrev[kind] = abbrev
            duplicates = Schema.get_duplicates(kind_to_abbrev.values())
        return kind_to_abbrev

    def create_abbreviations(self):
        """Creates abbreviations for node and edge kinds."""
        kind_to_abbrev = MetaGraph.find_abbrevs(self.node_dict.keys())
        kind_to_abbrev = {kind: abbrev.upper()
                          for kind, abbrev in kind_to_abbrev.items()}
        
        edge_set_to_keys = dict()
        for edge in self.edge_dict.keys():
            key = frozenset(map(str.lower, edge[:2]))
            value = edge[2]
            edge_set_to_keys.setdefault(key, list()).append(value)
        
        for edge_set, keys in edge_set_to_keys.items():
            key_to_abbrev = MetaGraph.find_abbrevs(keys)
            for key, abbrev in key_to_abbrev.items():
                previous_abbrev = kind_to_abbrev.get(key)
                if previous_abbrev and len(abbrev) <= len(previous_abbrev):
                    continue
                kind_to_abbrev[key] = abbrev
        
        self.set_abbreviations(kind_to_abbrev)
        self.kind_to_abbrev = kind_to_abbrev
        return kind_to_abbrev

    def set_abbreviations(self, kind_to_abbrev):
        for kind, node in self.node_dict.iteritems():
            node.abbrev = kind_to_abbrev[kind]
        for metaedge in self.edge_dict.itervalues():
            metaedge.kind_abbrev = kind_to_abbrev[metaedge.kind]
        
    def get_node(self, kind):
        return self.node_dict[kind]

    def get_edge(self, edge_tuple):
        return self.edge_dict[edge_tuple]

    def add_node(self, kind):
        metanode = MetaNode(kind)
        self.node_dict[kind] = metanode

    def add_edge(self, edge_tuple):
        """source_kind, target_kind, kind, direction"""
        source_kind, target_kind, kind, direction = edge_tuple
        source = self.get_node(source_kind)
        target = self.get_node(target_kind)
        
        metaedge = MetaEdge(source, target, kind, direction)
        self.edge_dict[edge_tuple] = metaedge
        source.edges.add(metaedge)

        if source != target or direction != 'both':
            inverse_direction = MetaEdge.inverse_direction[direction]
            inverse = MetaEdge(target, source, kind, inverse_direction)
            inverse_tuple = target_kind, source_kind, kind, inverse_direction
            self.edge_dict[inverse_tuple] = inverse
            target.edges.add(inverse)
            metaedge.inverse = inverse
            inverse.inverse = metaedge
        else:
            self.inverse = self
       
    def extract_metapaths(self, source_kind, target_kind, max_length):
        source = self.node_dict[source_kind]
        target = self.node_dict[target_kind]
        
        metapaths = ()
        metapaths = [MetaPath((edge,)) for edge in source.edges]
        previous_metapaths = list(metapaths)
        for depth in range(1, max_length):
            current_metapaths = list()
            for metapath in previous_metapaths:
                for add_edge in metapath.target().edges:
                    new_metapath = self.get_metapath(metapath.edges + (add_edge, ))
                    current_metapaths.append(new_metapath)
            metapaths.extend(current_metapaths)
            previous_metapaths = current_metapaths
        metapaths = [metapath for metapath in metapaths if metapath.target() == target]
        return metapaths
            
    def get_metapath(self, edges):
        """ """
        try:
            return self.path_dict[edges]
        except KeyError:
            assert isinstance(edges, tuple)
            metapath = MetaPath(edges)
            self.path_dict[edges] = metapath
            return metapath
        
class MetaNode(object):
    
    def __init__(self, kind):
        """ """
        self.kind = kind
        self.edges = set()
    
    def __hash__(self):
        return hash(self.kind)
    
    def __eq__(self, other):
        return self.kind == other.kind

    def __repr__(self):
        return self.kind

class MetaEdge(object):

    inverse_direction = {'forward': 'backward',
                         'backward': 'forward',
                         'both': 'both'}
    direction_to_abbrev = {'forward': '>', 'backward': '<', 'both': '-'}

    def __init__(self, source, target, kind, direction):
        """source and target are MetaNodes."""
        self.__dict__.update(locals())
        self.hash_ = hash((source, target, kind, direction))

    def __hash__(self):
        return self.hash_

    def __repr__(self):
        direction_abbrev = MetaEdge.direction_to_abbrev[self.direction]
        return '%s %s %s %s %s' % (self.source, direction_abbrev, self.kind, direction_abbrev, self.target)


class MetaPath(BasePath):
    
    def __init__(self, edges):
        """metaedges is a tuple of edges"""
        assert all(isinstance(edge, MetaEdge) for edge in edges)
        super(MetaPath, self).__init__(edges)

    def __repr__(self):
        s = ''
        for edge in self:
            source_abbrev = edge.source.abbrev
            dir_abbrev = MetaEdge.direction_to_abbrev[edge.direction]
            kind_abbrev = edge.kind_abbrev
            s += '{0}{1}{2}{1}'.format(source_abbrev, dir_abbrev, kind_abbrev)
        s+= self.target().abbrev
        return s


class Graph(object):
    
    def __init__(self, metagraph, data=dict()):
        """ """
        self.metagraph = metagraph
        self.data = data
        self.node_dict = dict()
        self.edge_dict = dict()

    def add_node(self, id_, kind, data=dict()):
        """ """
        metanode = self.metagraph.node_dict[kind]
        node = Node(id_, metanode, data)
        self.node_dict[id_] = node
        return node
    
    def add_edge(self, source_id, target_id, kind, direction, data=dict()):
        """ """
        source = self.node_dict[source_id]
        target = self.node_dict[target_id]
        metaedge_id = source.metanode.kind, target.metanode.kind, kind, direction
        metaedge = self.metagraph.edge_dict[metaedge_id]
        edge = Edge(source, target, metaedge, data)
        self.edge_dict[edge.get_id()] = edge
        
        inverse = Edge(target, source, metaedge.inverse, data)
        inverse_id = inverse.get_id()
        self.edge_dict[inverse_id] = inverse
        
        return edge, inverse

    def paths_from(self, source_id, metapath, no_duplicates=True):
        """ """
        source = self.node_dict[source_id]
        
        paths = list(Path((edge, )) for edge in source.edges[metapath[0]])

        for i in range(1, len(metapath)):
            current_paths = list()
            metaedge = metapath[i]
            for path in paths:
                nodes = path.get_nodes()
                edges = path.target().edges[metaedge]
                for edge in edges:
                    if no_duplicates and edge.target in nodes:
                        continue
                    newpath = path + Path((edge, ))
                    current_paths.append(newpath)
            paths = current_paths
        return paths
    
    def paths_between(self, source_id, target_id, metapath, no_duplicates=True):
        """can potentially divide paths"""
    
class Node(object):
    
    def __init__(self, id_, metanode, data):
        """ """
        self.__dict__.update(locals())
        self.edges = {metaedge: set() for metaedge in metanode.edges}

    def __hash__(self):
        return hash(self.id_)
    
    def __eq__(self, other):
        return self.id_ == other.id_

    def __repr__(self):
        return self.id_

class Edge(object):
    
    def __init__(self, source, target, metaedge, data):
        """source and target are Node objects. metaedge is the MetaEdge object
        representing the edge
        """
        assert isinstance(source, Node)
        assert isinstance(target, Node)
        assert isinstance(metaedge, MetaEdge)
        self.__dict__.update(locals())
        self.source.edges[metaedge].add(self)
        
    def get_id(self):
        return self.source.id_, self.target.id_, self.metaedge.kind, self.metaedge.direction

    def __hash__(self):
        return hash(self.get_id())
    
    def __eq__(self, other):
        return self.get_id() == other.get_id()
    
    def __repr__(self):
        source, target, kind, direction = self.get_id()
        direction = MetaEdge.direction_to_abbrev[direction]
        return '{0} {3} {2} {3} {1}'.format(source, target, kind, direction)

class Path(BasePath):
    
    def __init__(self, edges):
        """potentially metapath should be an input although it can be calculated"""
        super(Path, self).__init__(edges)
    
    def append(self, edge):
        self.edges.append(edge)
        self.nodes.append(edge.target)

    def __repr__(self):
        s = ''
        for edge in self:
            dir_abbrev = MetaEdge.direction_to_abbrev[edge.metaedge.direction]
            kind_abbrev = edge.metaedge.kind_abbrev
            s += '{0} {1} {2} {1} '.format(edge.source, dir_abbrev, edge.metaedge.kind)
        s = '{}{}'.format(s, self.target())
        return s

if __name__ == '__main__':
    """ """
    metaedges = [('gene', 'disease', 'association', 'both'),
             ('gene', 'gene', 'function', 'both'),
             ('gene', 'tissue', 'expression', 'both'),
             ('disease', 'tissue', 'pathology', 'both'),
             ('gene', 'gene', 'transcription', 'forward')]
    metagraph = MetaGraph.from_edge_tuples(metaedges)

    graph = Graph(metagraph)
    graph.add_node('IL17', 'gene')
    graph.add_node('BRCA1', 'gene')
    graph.add_node('ANAL1', 'gene')
    graph.add_node('MS', 'disease')
    graph.add_node('brain', 'tissue')
    
    graph.add_edge('MS', 'IL17', 'association', 'both')
    graph.add_edge('MS', 'brain', 'pathology', 'both')
    graph.add_edge('IL17', 'brain', 'expression', 'both')
    graph.add_edge('IL17', 'BRCA1', 'transcription', 'backward')

    
    metapaths = metagraph.extract_metapaths('gene', 'disease', 2)
    
    print graph.paths_from('BRCA1', metapaths[3])
    
    #print graph.node_dict.values()[0].edges
    #print graph.edge_dict

    #print nodes and edges
    #
    #print metapaths
