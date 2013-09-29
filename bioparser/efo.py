import os
import csv

import networkx

import data
import networkx_ontology

class EFO(networkx_ontology.NXObo):
    
    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('efo')
        obo_filename = 'efo.obo'
        json_filename = 'efo.networkx.json'
        pkl_filename = 'efo.networkx.pkl'
        self.keep_attributes = ['name', 'def_', 'synonym']
        super(EFO, self).__init__(directory, obo_filename, json_filename, pkl_filename)


    def get_graph_from_obo(self):
        """
        Store a parsed EFO obo file in self.ontology and a directed 
        networkx multigraph representation in self.graph.
        Colons are replaced with underscores node ids.
        """
        graph = super(EFO, self).get_graph_from_obo()
        #replace_colon = lambda s: s.replace(':', '_')
        #mapping = {node: replace_colon(node) for node in graph.nodes()}
        #networkx.relabel_nodes(graph, mapping, copy=False)
        return graph
    
    def get_diseases(self, root = 'EFO:0000408'):
        'EFO:0000408' # disease
        ontology = self.get_nx_ontology()
        return ontology.get_descendents(root)

    def get_neoplasms(self):
        root = 'EFO:0000616' # neoplasm
        ontology = self.get_nx_ontology()
        return ontology.get_descendents(root)
    
    def get_non_neoplastic_diseases(self):
        return self.get_diseases() - self.get_neoplasms()
    
    def gxa_query_compounds(self, root='CHEBI:37577'):
        """
        Finds all leaf nodes which are descendents of root.
        Colons are replaced with underscores in the returned output.
        Writes 'gxa_query_compounds.txt' showing term id and names.
        The default root is CHEBI:37577 for 'chemical compound'.
        """
        self.get_graph()        
        chemical_compounds = list(networkx.dfs_postorder_nodes(self.graph, source=root)) 
        query_compounds = filter(lambda x: self.graph.out_degree(x) == 0, chemical_compounds)
        self.write_terms(query_compounds, 'gxa_query_compounds.txt')
        #query_compounds = map(replace_colon, query_compounds)
        return query_compounds
    
    def gxa_query_diseases(self, root='EFO:0000408'):
        """
        The sef of leaf nodes which are descendents of root is computed.
        The predecessors of the leaf nodes are then computed. Nodes which have 
        more than five total descedents are excluded. Nodes that are descendents
        of included nodes are excluded.
        Writes 'gxa_query_diseases.txt' showing term id and names.
        The default root is EFO_0000408 for 'disease'.
        """
        self.get_graph()
        diseases = list(self.get_diseases(root))
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


if __name__ =='__main__':
    efo = EFO()
    efo.gxa_query_compounds()
    efo.gxa_query_diseases()