import re
import simplejson as json
import collections
import os
import gzip

import hetnet

def read_json(path):
    """ """
    
    if path.endswith('.gz'):
        read_file = gzip.open(path)
    else:
        read_file = open(path)
    json_object = json.load(read_file)
    read_file.close()
    
    metaedge_tuples = json_object['metaedge_tuples']
    metaedge_tuples = map(tuple, metaedge_tuples)
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)

    nodes = json_object['nodes']
    for node in nodes:
        graph.add_node(**node)

    edges = json_object['edges']
    for edge in edges:
        graph.add_edge(**edge)
    
    return graph

def write_json(graph, path):
    """ """
    metanode_kinds, metaedge_tuples, nodes, edges = get_graph_writables(graph)

    combined = collections.OrderedDict()
    combined['metanode_kinds'] = metanode_kinds
    combined['metaedge_tuples'] = metaedge_tuples
    combined['nodes'] = nodes
    combined['edges'] = edges


    if path.endswith('.gz'):
        write_file = gzip.open(path, 'w')
    else:
        write_file = open(path, 'w')
    json.dump(combined, write_file, indent=2)
    write_file.close()

def get_graph_writables(graph):
    """ """
    metanode_kinds = graph.metagraph.node_dict.keys()
    
    metaedge_tuples = [edge.get_id() for edge in
                       graph.metagraph.get_edges(exclude_inverts=True)]
    
    nodes = list()
    for node in graph.node_dict.itervalues():
        node_as_dict = collections.OrderedDict()
        node_as_dict['id_'] = node.id_
        node_as_dict['kind'] = node.metanode.id_
        node_as_dict['data'] = node.data
        nodes.append(node_as_dict)

    edges = list()
    for edge in graph.get_edges(exclude_inverts=True):
        edge_id_keys = ('source_id', 'target_id', 'kind', 'direction')
        edge_id = edge.get_id()
        edge_items = zip(edge_id_keys, edge_id)
        edge_as_dict = collections.OrderedDict(edge_items)
        edge_as_dict['data'] = edge.data
        edges.append(edge_as_dict)

    return metanode_kinds, metaedge_tuples, nodes, edges

def write_gml_graph(graph, path):
    """ """
    metanode_kinds, metaedge_tuples, nodes, edges = get_graph_writables(graph)
    
    with open(path, 'w') as write_file:
        gml_writer = GMLWriter(write_file)
        gml_writer.write_graph(nodes, edges)


class GMLWriter(object):
    """
    http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
    """
    
    def __init__(self, write_file):
        """GML writing and reading class"""
        self.gml_file = write_file # file to write GML to        
        self.write_indent = '\t'
        self.write_level = 0 # indentation level while writing
        
    def write_graph(self, nodes, edges):
        """nodes and edges are lists of dictionaries."""
        
        with GMLBlock(self, 'graph'):
            
            for node in nodes:
                with GMLBlock(self, 'node'):
                    self.write_properties(node)
                    
            for edge in edges:
                with GMLBlock(self, 'edge'):
                    self.write_properties(edge)
    

    def write(self, s):
        """Write string s to self.gml_file prepending the proper indentation."""
        indent = self.write_indent * self.write_level
        self.gml_file.write(indent + s)

    def write_properties(self, dictionary):
        for key, value in dictionary.items():
            self.write_property(key, value)

    def write_property(self, key, value, printing = False):
        """ """
        if not re.match(r'[A-Za-z]\w*\Z', key):
            if printing: print 'Invalid Key:', key
            return
        if isinstance(value, (int, long, float)):
            value = str(value)
        
        elif isinstance(value, basestring):
            #value = value.replace('"', "'")
            #value = value.replace('&', "AMPERSAND")
            if re.search(r'[&"\\]', value):
                if printing: print 'Invalid Value:', value
                return
        
        elif isinstance(value, (list, tuple, set)):
            with GMLBlock(self, key):
                for elem in value:
                    self.write_property('list', elem)
            return
        
        elif isinstance(value, dict):
            with GMLBlock(self, key):
                self.write_properties(value)
            return
        
        else:
            print 'GML formating not specified for', type(value)
            return

        line = '{} "{}"\n'.format(key, value)
        if len(line) > 254:
            if printing: print 'Line too long:', line
            return
        self.write(line)

class GMLBlock(object):
    
    def __init__(self, gml, key):
        self.gml = gml
        self.key = key
    
    def __enter__(self):
        self.gml.write('%s [\n' % self.key)
        self.gml.write_level += 1
    
    def __exit__(self, *args, **kwargs):
        self.gml.write_level -= 1
        self.gml.write(']\n')
