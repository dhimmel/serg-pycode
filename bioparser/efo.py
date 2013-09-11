import os
import csv
import json

import networkx
import networkx.readwrite.json_graph

import data
import obo


class EFO(object):
    
    def __init__(self, efo_dir=None):
        if efo_dir is None:
            efo_dir = data.current_path('efo')
        self.efo_dir = efo_dir
        self.obo_path = os.path.join(efo_dir, 'efo.obo')
        self.graph_json_path = os.path.join(efo_dir, 'efo.networkx.json')
        self.graph_pkl_path = os.path.join(efo_dir, 'efo.networkx.pkl')
        

    def get_graph_from_obo(self):
        """
        Store a parsed EFO obo file in self.ontology and a directed 
        networkx multigraph representation in self.graph.
        Colons are replaced with underscores node ids.
        """
        ontology = obo.OBOOntology(self.obo_path)
        graph = ontology.to_directed_networkx()
        replace_colon = lambda s: s.replace(':', '_')
        mapping = {node: replace_colon(node) for node in graph.nodes()}
        networkx.relabel_nodes(graph, mapping, copy=False)
        return graph

    def get_graph(self):
        """Returns the networkx graph representation of the EFO."""
        if not hasattr(self, 'graph'):
            if os.path.exists(self.graph_pkl_path):
                # Read networkx graph from json file
                with open(self.graph_json_path) as json_file:
                    serialized = json.load(json_file)
                # Networkx versions before 1.8 do not save edge keys in json serializations
                #self.graph = networkx.readwrite.json_graph.node_link_graph(serialized, True, True) 
                self.graph = networkx.read_gpickle(self.graph_pkl_path)

            else:
                # Write networkx graph to json file
                self.graph = self.get_graph_from_obo()
                serialized = networkx.readwrite.json_graph.node_link_data(self.graph)
                with open(self.graph_json_path, 'w') as json_file:
                    json.dump(serialized, json_file, indent=4)
                networkx.write_gpickle(self.graph, self.graph_pkl_path)

        return self.graph
    
    def get_id_to_name(self):
        graph = self.get_graph()
        id_to_name = dict()
        for node, data in graph.nodes_iter(data=True):
            id_to_name[node] = data['name']
        return id_to_name
    
    def get_descendents(self, node):
        return set(networkx.dfs_postorder_nodes(self.graph, source=node))
    
    def get_diseases(self):
        root = 'EFO_0000408' # disease
        return self.get_descendents(root)

    def get_neoplasms(self):
        root = 'EFO_0000616' # neoplasm
        return self.get_descendents(root)
    
    def get_non_neoplastic_diseases(self):
        return self.get_diseases() - self.get_neoplasms()
    
    def gxa_query_compounds(self, root='CHEBI_37577'):
        """
        Finds all leaf nodes which are descendents of root.
        Colons are replaced with underscores in the returned output.
        Writes 'gxa_query_compounds.txt' showing term id and names.
        The default root is CHEBI_37577 for 'chemical compound'.
        """
        self.get_graph()        
        chemical_compounds = list(networkx.dfs_postorder_nodes(self.graph, source=root)) 
        query_compounds = filter(lambda x: self.graph.out_degree(x) == 0, chemical_compounds)
        self.write_terms(query_compounds, 'gxa_query_compounds.txt')
        #query_compounds = map(replace_colon, query_compounds)
        return query_compounds
    
    def gxa_query_diseases(self, root='EFO_0000408'):
        """
        The sef of leaf nodes which are descendents of root is computed.
        The predecessors of the leaf nodes are then computed. Nodes which have 
        more than five total descedents are excluded. Nodes that are descendents
        of included nodes are excluded.
        Writes 'gxa_query_diseases.txt' showing term id and names.
        The default root is EFO_0000408 for 'disease'.
        """
        self.get_graph()
        diseases = list(self.get_disease_ids(root))
        leaf_diseases = filter(lambda x: self.graph.out_degree(x) == 0, diseases)
        query_diseases = set(leaf_diseases)
        for leaf_disease in leaf_diseases:
            query_diseases |= set(self.graph.predecessors(leaf_disease))
        
        for node in list(query_diseases):
            num_descendents = len(list(networkx.dfs_postorder_nodes(self.graph, source=node))) - 1
            if num_descendents > 5:
                query_diseases.remove(node)
        
        for node in list(query_diseases):
            descendents = set(networkx.dfs_postorder_nodes(self.graph, source=node))
            descendents.remove(node)
            query_diseases -= descendents
        
        self.write_terms(query_diseases, 'gxa_query_diseases.txt')
        #query_diseases = map(replace_colon, query_diseases)
        return query_diseases

    def write_terms(self, terms, file_name):
        path = os.path.join(self.efo_dir, file_name)
        f = open(path, 'w')
        writer = csv.writer(f, delimiter='\t')
        for term in terms:
            writer.writerow([term, self.graph.node[term]['name']])
        f.close()


if __name__ =='__main__':
    efo = EFO()
    efo.gxa_query_compounds()
    efo.gxa_query_diseases()