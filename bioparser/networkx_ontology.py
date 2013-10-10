import os
import csv
import json
import math

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
        return networkx.descendants(self.graph, node)
    
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

    def pairwise_shortest_paths(self, cutoff):
        node_set_to_length = dict()
        lengths = networkx.all_pairs_shortest_path_length(graph, cutoff)
        for source, lengths in lengths.iteritems():
            for target, length in lengths.iteritems():
                node_set = frozenset([source, target])
                node_set_to_length[fset] = length
        edges = list()
        for fset, length in node_set_to_length.iteritems():
            nodes = tuple(sorted(node_set))
            edge = nodes + (length, )
        return edges

    def annotate_information_content(self):
        """ """
        n_nodes = len(self.graph)
        for node, data in self.graph.nodes_iter(data=True):
            n_descendants = len(networkx.descendants(self.graph, node))
            node_ic = 1.0 - math.log1p(n_descendants) / math.log(n_nodes)
            data['IC'] = node_ic
    
    def intrinsic_similarity(self, node_0, node_1):
        """
        Returns a dictionary where the keys signify the method for 
        computing pairwise similarity and the values are the resulting
        similarities.
        
        Methods taken from "Seco N, Veale T, Hayes J. (2004) An Intrinsic
        Information Content Metric for Semantic Similarity in WordNet"
        """
        ic_0 = self.graph.node[node_0]['IC']
        ic_1 = self.graph.node[node_1]['IC']
        ancestors = (networkx.ancestors(self.graph, node_0) & 
                     networkx.ancestors(self.graph, node_1))
        resnik = max(self.graph.node[node]['IC'] for node in ancestors)
        lin = 2 * resnik / (ic_0 + ic_1)
        jen_distance = (ic_0 + ic_1) - 2 * resnik
        jen_similarity = 1 - jen_distance / 2
        metrics = {'resnik_similarity': resnik,
                   'lin_similarity': lin,
                   'jen_similarity': jen_similarity}
        return metrics
        

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
