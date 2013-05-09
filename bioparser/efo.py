import os
import csv

import networkx

import data
import obo

replace_colon = lambda s: s.replace(':', '_')

class EFO(object):
    
    def __init__(self, efo_dir=None):
        if efo_dir is None:
            efo_dir = data.current_path('efo')
        self.efo_dir = efo_dir
        self.obo_path = os.path.join(efo_dir, 'efo.obo.txt')

    def read(self):
        """
        Store a parsed EFO obo file in self.ontology and a directed 
        networkx multigraph representation in self.graph.
        """
        self.ontology = obo.OBOOntology(self.obo_path)
        self.graph = self.ontology.to_directed_networkx()

    def gxa_query_compounds(self, root='CHEBI:37577'):
        """
        Finds all leaf nodes which are descendents of root.
        Colons are replaced with underscores in the returned output.
        Writes 'gxa_query_compounds.txt' showing term id and names.
        The default root is CHEBI_37577 for 'chemical compound'.
        """
        if not hasattr(self, 'graph'):
            self.read()
        
        chemical_compounds = list(networkx.dfs_postorder_nodes(self.graph, source=root)) 
        query_compounds = filter(lambda x: self.graph.out_degree(x) == 0, chemical_compounds)
        self.write_terms(query_compounds, 'gxa_query_compounds.txt')
        query_compounds = map(replace_colon, query_compounds)
        return query_compounds

    def gxa_query_diseases(self, root='EFO:0000408'):
        """
        The sef of leaf nodes which are descendents of root is computed.
        The predecessors of the leaf nodes are then computed. Nodes which have 
        more than five total descedents are excluded. Nodes that are descendents
        of included nodes are excluded.
        Colons are replaced with underscores in the returned output.
        Writes 'gxa_query_diseases.txt' showing term id and names.
        The default root is EFO_0000408 for 'disease'.
        """
        if not hasattr(self, 'graph'):
            self.read()
        diseases = list(networkx.dfs_postorder_nodes(self.graph, source=root))
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
        query_diseases = map(replace_colon, query_diseases)
        return query_diseases

    def write_terms(self, terms, file_name):
        path = os.path.join(self.efo_dir, file_name)
        f = open(path, 'w')
        writer = csv.writer(f, delimiter='\t')
        for term in terms:
            writer.writerow([replace_colon(term), self.graph.node[term]['name']])
        f.close()


if __name__ =='__main__':
    efo = EFO()
    efo.gxa_query_compounds()
    efo.gxa_query_diseases()