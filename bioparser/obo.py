import json
import os
import re
import logging

import networkx

import networkx_ontology

singular_header_tags = ['format-version', 'data-version', 'date', 'saved-by']
singular_stanza_tags = ['id', 'name', 'is_anonymous', 'def', 'comment',
                 'is_obsolete', 'replaced_by']
plural_stanza_tags = ['alt_id', 'subset', 'synonym', 'xref', 'is_a']

trailing_modifiers_pattern = re.compile(r'\{(.*)\}$')

def parse_tag_value_pair(line):
    comment_split = line.split(' ! ', 1)
    tag_value_pair = comment_split[0]
    comment = comment_split[1] if len(comment_split) == 2 else None
    tag_value_pair = re.sub(trailing_modifiers_pattern, '', tag_value_pair)
    tag, value = tag_value_pair.split(': ', 1)
    return tag, value
    
def parse_obo(path):
    """
    http://www.geneontology.org/GO.format.obo-1_2.shtml
    """
    read_file = open(path)
    
    skip_line = lambda line: line.startswith('!') or not line
    
    # parse header section
    header = dict() # list of tag-value pairs
    line = read_file.next().rstrip('\r\n')
    while not line.startswith('['):
        line = line.rstrip('\r\n')
        if not skip_line(line):
            tag, value = parse_tag_value_pair(line)
            if tag in singular_header_tags:
                assert tag not in header
                header[tag] = value
            else:
                header.setdefault(tag, list()).append(value)

        line = next(read_file, None)

    # parse stanzas
    stanzas = list()
    while line is not None:
        line = line.rstrip('\r\n')
        stanza_name = line[1: -1]
        stanza = {'stanza_name': stanza_name}
        
        line = next(read_file, None)
        while line is not None and not line.startswith('['):
            line = line.rstrip('\r\n')
            if not skip_line(line):
                tag, value = parse_tag_value_pair(line)
                if tag in singular_stanza_tags:
                    if tag in stanza:
                        logging.warning('Multiple ' + tag + ' for ' + stanza['id'])
                    stanza[tag] = value
                else:
                    stanza.setdefault(tag, list()).append(value)
            line = next(read_file, None)
        stanzas.append(stanza)
    
    read_file.close()
    return header, stanzas

def process_stanzas(stanzas):
    for stanza in stanzas:
        if 'relationship' in stanza:
            relationship_dict = dict()
            for relationship in stanza['relationship']:
                type_id, term_id = relationship.split(' ', 1)
                relationship_dict.setdefault(type_id, list()).append(term_id)
            stanza['relationship'] = relationship_dict
    return stanzas

def terms_from_stanzas(stanzas):
    terms = list()
    for stanza in stanzas:
        if stanza['stanza_name'] != 'Term':
            continue
        if stanza.get('is_obsolete', False):
            continue
        terms.append(stanza)
    return terms
    

class OBO(object):
    
    def __init__(self, directory, obo_filename,
                 node_attributes=['name']):
        if obo_filename.endswith('.txt') or obo_filename.endswith('.obo'):
            name = obo_filename[: -4]
        else:
            name = obo_filename
        json_filename = '{}.json'.format(name)
        pkl_filename = '{}.networkx.pkl'.format(name)
        self.directory = directory
        self.obo_path = os.path.join(directory, obo_filename)
        self.graph_json_path = os.path.join(directory, json_filename)
        self.graph_pkl_path = os.path.join(directory, pkl_filename)
        self.node_attributes = node_attributes
    
    def to_networkx(self):
        """
        Creates a networkx directed multigraph representing the ontology.
        """
        header, stanzas = parse_obo(self.obo_path)
        process_stanzas(stanzas)
        terms = terms_from_stanzas(stanzas)
        
        # create graph
        graph = networkx.MultiDiGraph()
        graph.graph.update(header)
        
        # add nodes
        for stanza in terms:
            
            node = stanza['id']
            node_data = {key: stanza[key] for key in self.node_attributes if key in stanza}
            graph.add_node(node, node_data)
        
        # add edges
        for stanza in terms:
            child = stanza['id']
            type_to_terms = dict()
            if 'is_a' in stanza:
                type_to_terms['is_a'] = stanza['is_a']
            if 'relationship' in stanza:
                type_to_terms.update(stanza['relationship'])
            for key, parents in type_to_terms.items():
                for parent in parents:
                    graph.add_edge(parent, child, key=key)
                    
        
        
        assert networkx.is_directed_acyclic_graph(graph)
        return graph
    
    def write_graph(self):
        serialized = networkx.readwrite.json_graph.node_link_data(self.graph)
        with open(self.graph_json_path, 'w') as json_file:
            json.dump(serialized, json_file, indent=2)
        networkx.write_gpickle(self.graph, self.graph_pkl_path)
    
    def read_graph(self):
        return networkx.read_gpickle(self.graph_pkl_path)
    
    def get_graph(self):
        """Returns the networkx graph representation of the EFO."""
        if not hasattr(self, 'graph'):
            if os.path.exists(self.graph_pkl_path):
                self.graph = self.read_graph()
            else:
                self.graph = self.to_networkx()
                self.write_graph()

        return self.graph

    def get_ontology(self):
        """ """
        if not hasattr(self, 'ontology'):
            graph = self.get_graph()
            self.ontology = networkx_ontology.NXOntology(graph)
        return self.ontology



#stanza_name = re.search(r'^\[(.*)\]$', line).group(1)
if __name__ =='__main__':
    path = '/home/dhimmels/Documents/serg/data-sources/bto/131007/BrendaTissueOBO.txt'
    graph = to_networkx(path)
    print graph.nodes()