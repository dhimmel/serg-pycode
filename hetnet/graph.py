# daniel.himmelstein@gmail.com
import pickle as pickle

import readwrite

direction_to_inverse = {'forward': 'backward',
                         'backward': 'forward',
                         'both': 'both'}

direction_to_abbrev = {'forward': '>', 'backward': '<', 'both': '-'}


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
    
    def inverse_edges(self):
        return tuple(reversed(list(edge.inverse for edge in self)))
    
    def traverses_mask(self):
        """Whether any nodes or edges contained in the path are masked."""
        for edge in self:
            if edge.masked or edge.source.masked:
                return True
        return self.target().masked
    
    def __iter__(self):
        return iter(self.edges)

    def __getitem__(self, key):
        return self.edges[key]

    def __len__(self):
        return len(self.edges)

    def __hash__(self):
        return hash(self.edges)
    
    def __eq__(self):
        return self.edges == other.edges


class MetaGraph(object):
    
    def __init__(self):
        """ """
        self.node_dict = dict()
        self.edge_dict = dict()
        self.path_dict = dict()
    
    def __getstate__(self):
        state = self.__dict__.copy()
        state['node_dict_items'] = state.pop('node_dict').items()
        state['path_dict_items'] = state.pop('path_dict').items()
        state['edge_dict_items'] = state.pop('edge_dict').items()
        return state

    def __setstate__(self, state):
        state['node_dict'] = dict(state['node_dict_items'])
        state['edge_dict'] = dict(state['edge_dict_items'])
        state['path_dict'] = dict(state['path_dict_items'])
        self.__dict__.update(state)

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

        if source == target and direction == 'both':
            metaedge.inverse = metaedge
        else:
            inverse_direction = direction_to_inverse[direction]
            inverse = MetaEdge(target, source, kind, inverse_direction)
            inverse_tuple = target_kind, source_kind, kind, inverse_direction
            self.edge_dict[inverse_tuple] = inverse
            target.edges.add(inverse)
            metaedge.inverse = inverse
            inverse.inverse = metaedge

       
    def extract_metapaths(self, source_kind, target_kind, max_length):
        source = self.node_dict[source_kind]
        target = self.node_dict[target_kind]
        
        metapaths = [self.get_metapath((edge, )) for edge in source.edges]
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
            """
            inverse_edges = metapath.inverse_edges()
            inverse = MetaPath(inverse_edges)
            self.path_dict[inverse_edges] = inverse

            metapath.inverse = inverse
            inverse.inverse = metapath
            """
            return metapath
    
    
class MetaNode(object):
    
    def __init__(self, kind):
        """ """
        self.kind = kind
        self.edges = set()
        self.masked = False
    
    def __hash__(self):
        return hash(self.kind)

    def __getstate__(self):
        state = self.__dict__.copy()
        state['edges_list'] = list(state.pop('edges'))
        return state
    
    def __setstate__(self, state):
        state['edges'] = set(state['edges_list'])
        self.__dict__.update(state)    

    def __eq__(self, other):
        return self.kind == other.kind

    def __repr__(self):
        return self.kind

class MetaEdge(object):
    """ """
    

    def __init__(self, source, target, kind, direction):
        """source and target are MetaNodes."""
        #self.__dict__.update(locals())
        self.source = source
        self.target = target
        self.kind = kind
        self.direction = direction
        self.hash_ = hash(self.get_id())
        self.masked = False


    def __getstate__(self):
        state = self.__dict__.copy()
        state['inverse_id'] = state.pop('inverse').get_id()
        return state


    def get_id(self):
        """ """
        return self.source.kind, self.target.kind, self.kind, self.direction

    def __hash__(self):
        return self.hash_

    """
    def __eq__(self, other):
        self.get_id() == other.get_id()
    """
    
    def __repr__(self):
        dir_abbrev = direction_to_abbrev[self.direction]
        return '{0} {3} {2} {3} {1}'.format(
            self.source, self.target, self.kind, dir_abbrev)


class MetaPath(BasePath):
    
    def __init__(self, edges):
        """metaedges is a tuple of edges"""
        assert all(isinstance(edge, MetaEdge) for edge in edges)
        super(MetaPath, self).__init__(edges)
    
    def __repr__(self):
        s = ''
        for edge in self:
            source_abbrev = edge.source.abbrev
            dir_abbrev = direction_to_abbrev[edge.direction]
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
        edge.inverted = False
        
        inverse = Edge(target, source, metaedge.inverse, data)
        inverse_id = inverse.get_id()
        self.edge_dict[inverse_id] = inverse
        inverse.inverted = True
       
        return edge, inverse

    def paths_from(self, source, metapath,
                   duplicates=False, masked=True,
                   exclude_nodes=set(), exclude_edges=set()):
        """
        Return a list of Paths starting with source and following metapath.
        Setting duplicates False disallows paths with repeated nodes.
        Setting masked False disallows paths which traverse a masked node or edge.
        exclude_nodes and exclude_edges allow specification of additional nodes
        and edges beyond (or independent of) masked nodes and edges.
        """

        if not isinstance(source, Node):
            source = self.node_dict[source]
        
        if masked and source.masked:
            return None
        
        if source in exclude_nodes:
            return None
        
        paths = list()
        print metapath[0]
        print source.edges
        for edge in source.edges[metapath[0]]:
            edge_target = edge.target
            if edge_target in exclude_nodes:
                continue
            if edge in exclude_edges:
                continue
            if not masked and (edge_target.masked or edge.masked):
                continue
            if not duplicates and edge_target == source:
                continue
            path = Path((edge, ))
            paths.append(path)
        
        for i in range(1, len(metapath)):
            current_paths = list()
            metaedge = metapath[i]
            for path in paths:
                nodes = path.get_nodes()
                edges = path.target().edges[metaedge]
                for edge in edges:
                    edge_target = edge.target
                    if edge_target in exclude_nodes:
                        continue
                    if edge in exclude_edges:
                        continue
                    if not masked and (edge_target.masked or edge.masked):
                        continue
                    if not duplicates and edge_target in nodes:
                        continue
                    newpath = Path(path.edges + (edge, ))
                    current_paths.append(newpath)
            paths = current_paths
        
        return paths
    
    
    def paths_between(self, source, target, metapath,
                      duplicates=False, masked=True,
                      exclude_nodes=set(), exclude_edges=set()):
        """
        Retreive the paths starting with the node source and ending on the
        node target. Future implementations should split the metapath, computing
        paths_from the source and target and look for the intersection at the
        intermediary Node position.
        """
        paths = self.paths_from(source, metapath, duplicates, masked, exclude_nodes, excluded_edges)
        paths = [path for path in paths if path.target == target]
        return paths        
    
    def unmask(self):
        """Unmask all nodes and edges contained within the graph"""
        for dictionary in self.node_dict, self.edge_dict:
            for value in dictionary.itervalues():
                value.masked = False
    
    def write_gml(self, path):
        writer = readwrite.HetnetGMLWriter(self, path)
        writer.export_graph()

    #"""
    def __getstate__(self):
        state = self.__dict__.copy()
        state['node_dict_items'] = state.pop('node_dict').items()
        state['edge_dict_items'] = state.pop('edge_dict').items()
        return state

    def __setstate__(self, state):
        state['node_dict'] = dict(state['node_dict_items'])
        state['edge_dict'] = dict(state['edge_dict_items'])
        self.__dict__.update(state)
    #"""
    def write_pickle(self, path):
        """
        http://bugs.python.org/issue1761028
        http://bugs.python.org/issue1062277
        http://bugs.python.org/issue998998
        http://bugs.python.org/issue1581183
        http://bugs.python.org/issue9269
        """
        with open(path, 'w') as write_file:
            pickle.dump(self, write_file, protocol=-1)

    @staticmethod
    def read_pickle(path):
        """
        http://bugs.python.org/issue1761028
        http://bugs.python.org/issue1062277
        http://bugs.python.org/issue998998
        http://bugs.python.org/issue1581183
        """
        with open(path) as read_file:
            graph = pickle.load(read_file)


        # Reconstruct metagraph inverse references
        metaedge_dict = graph.metagraph.edge_dict
        for metaedge in metaedge_dict.itervalues():
            inverse_id = metaedge.__dict__.pop('inverse_id')
            inverse_metaedge = metaedge_dict[inverse_id]
            metaedge.inverse = inverse_metaedge
            print metaedge

        """
        # Reconstruct edge source and target Node references
        node_dict = graph.node_dict
        for edge in graph.edge_dict.itervalues():
            edge.source = node_dict[edge.source_id]
            edge.target = node_dict[edge.target_id]
        """
        return graph
        
    
class Node(object):
    
    def __init__(self, id_, metanode, data):
        """ """
        self.__dict__.update(locals())
        self.edges = {metaedge: set() for metaedge in metanode.edges}
        self.masked = False

    def get_edges(self, metaedge, exclude_masked=True):
        
        if exclude_masked:
            edges = list()
            for edge in self.edges[metaedge]:
                if edge.masked or edge.target.masked:
                    continue
                edges.append(edge)
        else:
            edges = self.edges[metaedge]
        return edges

    def __getstate__(self):
        state = self.__dict__.copy()
        state['edges_items'] = [(metaedge, list(edge_set))
            for metaedge, edge_set in state.pop('edges').items()]
        return state

    def __setstate__(self, state):
        state['edges'] = {metaedge: set(edge_list) for metaedge, edge_list in state['edges_items']}
        self.__dict__.update(state)


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
        self.hash_ = hash(self.get_id())
        self.source.edges[metaedge].add(self)
        self.masked = False
    
    def mask(self):
        self.masked = True
        
    def unmask(self):
        self.masked = False
        
    def get_id(self):
        return self.source.id_, self.target.id_, self.metaedge.kind, self.metaedge.direction

    """
    def __getstate__(self):
        state = self.__dict__
        state['source_id'] = state.pop('source').id_
        state['target_id'] = state.pop('target').id_
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
    """
    def __hash__(self):
        return self.hash_
    
    def __eq__(self, other):
        return self.get_id() == other.get_id()
    
    def __repr__(self):
        source, target, kind, direction = self.get_id()
        dir_abbrev = direction_to_abbrev[direction]
        return '{0} {3} {2} {3} {1}'.format(source, target, kind, dir_abbrev)

class Path(BasePath):
    
    def __init__(self, edges):
        """potentially metapath should be an input although it can be calculated"""
        super(Path, self).__init__(edges)
    
    def __repr__(self):
        s = ''
        for edge in self:
            dir_abbrev = direction_to_abbrev[edge.metaedge.direction]
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

    #metapaths = metagraph.extract_metapaths('gene', 'disease', 2)

    
    pkl_path = '/home/dhimmels/Desktop/test-graph.pkl'
    graph.write_pickle(pkl_path)
    graph = Graph.read_pickle(pkl_path)


    metapaths = metagraph.extract_metapaths('gene', 'disease', 2)
    print graph.paths_from('BRCA1', metapaths[3])

    """
    for node in graph.node_dict.values():
        print node.edges

    metapaths = metagraph.extract_metapaths('gene', 'disease', 2)
    """
    #print graph.node_dict.values()[0].edges
    #print graph.edge_dict

    #print nodes and edges
    #
    #print metapaths
