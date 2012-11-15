import collections
import types
import os
import re

# set defualt arguments for neo4j. Must execute before neo4j import.
os.environ['NEO4J_PYTHON_JVMARGS'] = '-Xms1G -Xmx2G -XX:+UseConcMarkSweepGC'
import neo4j

class Graph(object):
    
    node_kinds = [] # Must be set by user or derived class
    rel_kinds = [] # Must be set by user or derived class
    
    def __init__(self, neo_db_dir):
        """neo_db_dir contains the directory storing the neo4j database"""
        
        print ' Openning Database Connection '.center(79, '#')
        self.db = neo4j.GraphDatabase(neo_db_dir)
        
        if not self.db.node.indexes.exists('global'):
            
            print 'No global node index.'
            with self.db.transaction:
                self.db.node.indexes.create('global')
                self.db.relationship.indexes.create('global')              
            print 'global node index created'
            print 'global relationship index created'
    
    def shutdown(self):
        self.db.shutdown()
        print ' db connection shutdown '.center(79, '#')
        print ''.center(79, '#')

    def string_apply(self, elem, fxn):
        """Apply function fxn to strings"""
        if isinstance(elem, basestring):
            return fxn(elem)
        if isinstance(elem, types.GeneratorType):
            elem = list(elem)
        if isinstance(elem, list):
            return list(self.string_apply(e, fxn) for e in elem)
        if isinstance(elem, set):
            return set(self.string_apply(e, fxn) for e in elem)
        if isinstance(elem, tuple):
            return tuple(self.string_apply(e, fxn) for e in elem)
        if isinstance(elem, dict):
            return {self.string_apply(k, fxn): self.string_apply(v, fxn) for k, v in elem.items()}
        if not isinstance(elem, (int, long, float, bool, types.NoneType)):
            print type(elem)
        return elem

    def create_node(self, node_kind, properties, index=[]):
        """Create and return a node of kind node_kind in graph g. 
        properties contains the key, value pairs
        to create node properties. index is a list of the keys in properties
        that the user wants indexed in the node index 'global'.
        """
        
        string_fxn = lambda s: s.encode('ascii', 'replace')
        properties = self.string_apply(properties, string_fxn)
        
        with self.db.transaction:
            
            node = self.db.node(**properties)
            node['kind'] = node_kind
    
            node_idx = self.db.node.indexes.get('global')
            node_idx['kind'][node_kind] = node
                        
            for key in index:
                if key in properties:
                    value = properties[key]
                    if isinstance(value, basestring):
                        value = [value]
                    for element in value:
                        node_idx[key][element] = node
                
        return node

    def create_relationship(self, start, end, kind, properties=dict(), index=[]):
        """Create and return a relationship."""
        
        if self.get_relationship_between_nodes(start, end, kind):
            print kind, 'duplicate relationship'
            return

        with self.db.transaction:
            
            rel = start.relationships.create(kind, end, **properties)
            rel['kind'] = kind
            
            rel_idx = self.db.relationship.indexes.get('global')
            rel_idx['kind'][kind] = rel
            for key in index:
                if key in properties:
                    rel_idx[key][properties[key]] = rel
        
        return rel
    
    def get_relationship_between_nodes(self, start, end, kind):
        """Returns a relationship of a specific kind between two nodes if
        it exists. Direction dependent."""
        for rel in getattr(start, kind).outgoing:
            if rel.end == end:
                return rel
                
    
    def get_nodes_of_kind(self, kind):
        """Returns a list of nodes that are indexed having specified kind."""
        return self.query_nodes('kind:' + kind)
    
    def get_relationships_of_kind(self, kind):
        """Return relationships of the specified kind."""
        rel_idx = self.db.relationship.indexes.get('global')
        hits = rel_idx['kind'][kind]
        return list(hits)
    
    def query_nodes(self, query, index_name='global'):
        """Perform a lucene based query using a specified index to return
        a list of nodes meeting the query conditions.
        """
        try:
            query = query.encode('ascii', 'replace')
        except UnicodeDecodeError:
            print 'In QN: UnicodeDecodeError:', query
            return list()
        
        node_idx = self.db.node.indexes.get(index_name)
        hits = node_idx.query(query)
        return list(hits)

    
    def query_node(self, query, index_name='global', printing=False):
        """Meant to perform a lucene query when a single node is expected as
        the result. Either returns a single node or None if the query returns
        no results. If true, printing instructs the program to print queries
        that do not return exactly one hit.
        """
        node = None
        hits = self.query_nodes(query, index_name)
        num_of_hits = len(hits)
        if num_of_hits == 1:
            node = hits.pop()
        if printing:
            if num_of_hits > 1:
                print query, 'returns', num_of_hits, 'hits'
            elif num_of_hits < 1:
                print query, 'no hits'
        return node

    def print_node_counts(self):
        """Prints the number of each node kind in the graph."""
        print ' Node Type Counts '.center(79, '#')
        for kind in self.node_kinds:
            nodes = self.get_nodes_of_kind(kind)
            print len(nodes), kind, 'nodes'
    
    def print_relationship_counts(self):
        """Print the number of each relationship kind in the network."""
        print ' Relationship Type Counts '.center(79, '#')
        rel_kinds = [str(kind.toString()) for kind in self.db.relationshipTypes]
        for kind in rel_kinds:
            rels = self.get_relationships_of_kind(kind)
            print len(rels), kind, 'relationships'
    
    def delete_all_nodes_of_kind(self, kind):
        """Delete all instance nodes of the kind."""
        nodes = self.get_nodes_of_kind(kind)
        num_of_nodes = len(nodes)
        for node in nodes:
            self.delete_node(node)
        print 'All', num_of_nodes, kind, 'nodes deleted.'
    
    def delete_node(self, node):
        """Delete all relationships to or from the node and delete the node."""
        
        indexes = self.db.node.indexes
        idx_names = list(indexes._index.nodeIndexNames())
        
        with self.db.transaction:
            
            # Delete index references
            for idx_name in idx_names:
                idx = indexes.get(idx_name)
                del idx[node]
            
            # Delete relationships involving the node
            for rel in list(node.relationships):
                self.delete_relationship(rel)
           
            # Delete the node
            node.delete()
    
    def delete_relationship(self, rel):
        """Delete a single relationship and its index references"""
        
        indexes = self.db.relationship.indexes
        idx_names = list(indexes._index.relationshipIndexNames())
        
        with self.db.transaction:
            
            # Delete index references
            for idx_name in idx_names:
                idx = indexes.get(idx_name)
                del idx[rel]
            
            # Delete relationship
            rel.delete()
    
    def delete_all_relationships_of_kind(self, kind):
        """Delete relatioships of specified kind. Uses global relationship index"""
        relationships = self.get_relationships_of_kind(kind)
        rel_num = len(relationships)
        for rel in list(relationships):
            self.delete_relationship(rel)
        print 'All', rel_num, kind, 'relationships deleted.'
    
    def get_connected_nodes(self, node_kinds=None):
        if node_kinds is None:
            node_kinds = self.node_kinds
        if isinstance(node_kinds, basestring):
            raise Exception('node_kinds should be a non string iterable')
        nodes = list()
        for kind in node_kinds:
            nodes.extend(self.get_nodes_of_kind(kind))
        connected_nodes = filter(lambda node: len(node.relationships), nodes)
        print len(connected_nodes), 'of', len(nodes), 'nodes have a relationship'
        return connected_nodes
    
    def export_GML(self, nodes, rels, path):
        """
        Currently must use the following gremlin query
        g = new Neo4jGraph('/scratch/omicnet.db')
        g.saveGML('/home/dhimmels/Documents/serg/omicnet/omicnet.gml')
        g.shutdown()
        http://www.fim.uni-passau.de/fileadmin/files/lehrstuhl/brandenburg/projekte/gml/gml-technical-report.pdf
        """
        
        print 'exporting GML'
        level = 0
        gml = open(path, 'w')
        
        gml.write('graph' + ' [\n')
        for node in nodes:
            self.gml_write_elem(node, 'node', gml)
        for rel in rels:
            self.gml_write_elem(rel, 'edge', gml)

        gml.write(']')
        gml.close()
    
    def gml_write_elem(self, elem, kind, gml):
        """ """
        if kind in ['node', 'edge']:
            gml.write(kind + ' [\n')
            if kind == 'node':
                gml.write('id ' + str(elem.id) + '\n')
            elif kind == 'edge':
                gml.write('source ' + str(elem.start.id) + '\n')
                gml.write('target ' + str(elem.end.id) + '\n')
            for item in elem.items():
                self.gml_write_elem(item, 'other', gml)
            gml.write(']')
        # not a node or edge
        else:
            key, value = elem
            if not re.match(r'[A-Za-z]\w*\Z', key):
                print key
            if key == 'indication':
                value = 'fake indication'
            gml.write(key)
            if isinstance(value, (int, long, float)):
                gml.write(' ' + str(value))
            elif isinstance(value, basestring):
                gml.write(' "' + value + '"')
            else:
                gml.write(' "fake value"')
                """
                gml.write(': [\n')
                if isinstance(value, (list, tuple, set)):
                    for val in value:
                        gml.write('x "' + val + '"\n')
                if isinstance(value, dict):
                    for k, val in value:
                        gml.write(k + ' "' + val + '"\n')
                gml.write(']')
                """
        gml.write('\n')
            
        
            
    def create_subgraph(self, node_kinds, rel_kinds, sub_g, connected_nodes_only=True):
        """INCOMPLETE"""
        rels = list()
        for kind in rel_kinds:
            rel_of_kind = self.get_relationships_of_kind(kind)
            for rel in rel_of_kind:
                rels.extend()
        rels = set(rels)
        
        nodes = list()
        for kind in node_kinds:
            nodes_of_kind = self.get_nodes_of_kind(kind)
            nodes = filter(lambda n: n, nodes)
            nodes.extend(self.get_nodes_of_kind(kind))
        nodes = set(nodes)
        
        
        sub_g.shutdown()
        return sub_g
        
        