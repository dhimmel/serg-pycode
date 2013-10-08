import os
import csv
import json

import networkx
import networkx.readwrite.json_graph

import obo

class NXOntology(object):
    
    def __init__(self, graph):
        assert networkx.is_directed_acyclic_graph(graph)
        self.graph = graph

    def get_id_to_name(self):
        id_to_name = dict()
        for node, data in self.graph.nodes_iter(data=True):
            id_to_name[node] = data['name']
        return id_to_name
    
    def get_descendents(self, node):
        return set(networkx.dfs_postorder_nodes(self.graph, source=node))
    
    def initialize_attribute(self, attribute):
        for node, data in self.graph.nodes_iter(data=True):
            data[attribute] = set()

    def propogate_annotation(self, node, attribute, element):
        """
        Add element to the node attribute specified by attribute for 
        the specified node and all of its descencdants. The node attribute
        is assumed to be a set.
        """
        descendents = self.get_descendents(node)
        descendents.add(node)
        for descendent in descendents:
            data = self.graph.node[descendent]
            data.setdefault(attribute, set()).add(element)


class NXObo(object):
    
    def __init__(self, directory, obo_filename, json_filename, pkl_filename):
        self.directory = directory
        self.obo_path = os.path.join(directory, obo_filename)
        self.graph_json_path = os.path.join(directory, json_filename)
        self.graph_pkl_path = os.path.join(directory, pkl_filename)

    def get_graph_from_obo(self):
        """
        Store a parsed EFO obo file in self.ontology and a directed 
        networkx multigraph representation in self.graph.
        """
        ontology = obo.OBOOntology(self.obo_path)
        graph = ontology.to_directed_networkx(attributes=self.keep_attributes)
        return graph

    def get_graph(self):
        """Returns the networkx graph representation of the EFO."""
        if not hasattr(self, 'graph'):
            if os.path.exists(self.graph_pkl_path):
                # Read networkx graph from json file
                #with open(self.graph_json_path) as json_file:
                #    serialized = json.load(json_file)
                # Networkx versions before 1.8 do not save edge keys in json serializations
                #self.graph = networkx.readwrite.json_graph.node_link_graph(serialized, True, True) 
                self.graph = networkx.read_gpickle(self.graph_pkl_path)

            else:
                # Write networkx graph to json file
                self.graph = self.get_graph_from_obo()
                serialized = networkx.readwrite.json_graph.node_link_data(self.graph)
                with open(self.graph_json_path, 'w') as json_file:
                    json.dump(serialized, json_file, indent=2)
                networkx.write_gpickle(self.graph, self.graph_pkl_path)

        return self.graph

    def get_nx_ontology(self):
        if not hasattr(self, 'nx_ontology'):
            graph = self.get_graph()
            self.nx_ontology = NXOntology(graph)
        return self.nx_ontology
    
    def write_terms(self, terms, file_name):
        path = os.path.join(self.directory, file_name)
        f = open(path, 'w')
        writer = csv.writer(f, delimiter='\t')
        for term in terms:
            writer.writerow([term, self.graph.node[term]['name']])
        f.close()
