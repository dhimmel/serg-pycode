import os
import csv
import json
import math
import itertools

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
        the specified node and all of its descendants. The node attribute
        is assumed to be a set.
        """
        descendents = self.get_descendents(node)
        descendents.add(node)
        for descendent in descendents:
            data = self.graph.node[descendent]
            data.setdefault(attribute, set()).add(element)

    def get_roots(self):
        roots = set()
        for node in self.graph.nodes_iter():
            if not(self.graph.in_edges(node)):
                roots.add(node)
        return roots

    def most_specific_superior_annotations(self, attribute, new_attribute):
        """
        For each paths from a node to a root, find the most specific annotated node and consider
        those annotations as members of the most_specific_superior_annotations new_attribute.
        """
        self.initialize_attribute(new_attribute)
        memo_dict = dict()
        roots = self.get_roots()

        for node, data in self.graph.nodes_iter(data=True):
            superior_annotation = set()
            for root in roots:
                paths = self.all_paths(node, root, memo_dict)
                for path in paths:
                    direct_annotation = set()
                    for elem in reversed(path):
                        direct_annotation = self.graph.node[elem][attribute]
                        if direct_annotation:
                            break
                    superior_annotation.update(direct_annotation)
            data[new_attribute] = superior_annotation


    def all_paths(self, inferior_node, superior_node, memo_dict = None):
        if memo_dict is None:
            # putting {}, or any other mutable object
            # as the default argument is wrong
            memo_dict = dict()

        if inferior_node == superior_node: # Don't memoize trivial case
            return frozenset([(inferior_node,)])
        else:
            pair = (inferior_node, superior_node)
            if pair in memo_dict: # Is answer memoized already?
                return memo_dict[pair]
            else:
                result = set()
                for new_source in self.graph.predecessors(inferior_node):
                    paths = self.all_paths(new_source, superior_node, memo_dict)
                    for path in paths:
                        path = (inferior_node,) + path
                        result.add(path)
                result = frozenset(result)
                # Memoize answer
                memo_dict[(inferior_node, superior_node)] = result
                return result


    def pairwise_shortest_paths(self, cutoff):
        """
        Calculates node pairs with a shortest path distance less than or equal to cutoff.
        A list of tuples is returned where each tuple represents a node pair
        and is formatted (node_1, node_2, distance).
        """
        node_set_to_length = dict()
        lengths = networkx.all_pairs_shortest_path_length(self.graph, cutoff)
        for source, lengths in lengths.iteritems():
            for target, length in lengths.iteritems():
                node_set = frozenset([source, target])
                node_set_to_length[fset] = length
        edges = list()
        for fset, length in node_set_to_length.iteritems():
            nodes = tuple(sorted(node_set))
            edge = nodes + (length, )
        return edges

    def annotate_instrinsic_IC(self):
        """ """
        n_nodes = len(self.graph)
        for node, data in self.graph.nodes_iter(data=True):
            n_descendants = len(networkx.descendants(self.graph, node))
            node_ic = 1.0 - math.log1p(n_descendants) / math.log(n_nodes)
            data['intrinsic_IC'] = node_ic
    
    def intrinsic_similarity(self, node_0, node_1):
        """
        Returns a dictionary where the keys signify the method for 
        computing pairwise similarity and the values are the resulting
        similarities.
        
        Methods taken from "Seco N, Veale T, Hayes J. (2004) An Intrinsic
        Information Content Metric for Semantic Similarity in WordNet"
        """
        ic_0 = self.graph.node[node_0]['intrinsic_IC']
        ic_1 = self.graph.node[node_1]['intrinsic_IC']
        subsume_0 = networkx.ancestors(self.graph, node_0)
        subsume_1 = networkx.ancestors(self.graph, node_1)
        subsume_0.add(node_0)
        subsume_1.add(node_1)
        subsume = subsume_0 & subsume_1
        resnik = max(self.graph.node[node]['intrinsic_IC'] for node in subsume)
        lin = 2 * resnik / (ic_0 + ic_1)
        jen_distance = (ic_0 + ic_1) - 2 * resnik
        jen_similarity = 1 - jen_distance / 2
        metrics = {'resnik_similarity': resnik,
                   'lin_similarity': lin,
                   'jen_similarity': jen_similarity}
        return metrics

    def pairwise_similarities(self, nodes):
        self.annotate_instrinsic_IC()
        for node_0, node_1 in itertools.combinations(nodes, 2):
            metrics = self.intrinsic_similarity(node_0, node_1)
            metrics['source'] = node_0
            metrics['target'] = node_1
            yield metrics

