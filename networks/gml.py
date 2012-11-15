import collections
import re

class Writer(object):
    """
    Alternative gremlin query
    g = new Neo4jGraph('/scratch/omicnet.db')
    g.saveGML('/home/dhimmels/Documents/serg/omicnet/omicnet.gml')
    g.shutdown()

    http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
    """
    
    def __init__(self, path):
        """GML writing and reading class"""
        self.gml_file = open(path, 'w') # file to write GML to        
        self.write_indent = '\t'
        self.write_level = 0 # indentation level while writing
        
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
            value = '"%s"' % value
        
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
            print type(value)
            return
        
        line = '%s %s\n' % (key, value)
        if len(line) > 254:
            if printing: print 'Line too long:', line
            return
        self.write(line)

class Neo4jWriter(Writer):
    
    def __init__(self, path):
        super(Neo4jWriter, self).__init__(path)
        #self.edge_num = 0
    
    def export_graph(self, nodes, rels):
        """Write gml file from a neo4j graph set of nodes and rels."""
        with GMLBlock(self, 'graph'):
            for node in nodes:
                self.write_node(node)
            for rel in rels:
                self.write_edge(rel)
        self.gml_file.close()

    def write_node(self, node):
        """ """
        with GMLBlock(self, 'node'):
            self.write_property('id', node.id)
            self.write_properties(node)
    
    def write_edge(self, rel):
        """ """
        with GMLBlock(self, 'edge'):
            #self.write_property('id', self.edge_num)
            #self.edge_num += 1
            self.write_property('source', rel.start.id)
            self.write_property('target', rel.end.id)
            self.write_properties(rel)
        
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

