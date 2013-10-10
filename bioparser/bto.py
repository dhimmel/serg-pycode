import os
import csv

import networkx

import data
import obo
import networkx_ontology

class BTO(obo.OBO):
    
    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('bto')
        obo_filename = 'BrendaTissueOBO.txt'
        node_attributes = ['name', 'def_', 'synonym']
        relationship_tags = ['is_a', 'part_of', 'develops_from']
        super(BTO, self).__init__(directory, obo_filename, node_attributes)

    def get_animal_tissues(self):
        root = 'BTO:0001489' # whole body (descendant of animal)
        nodes = networkx.descendants(self.get_graph(), root)
        nodes.add(root)
        return nodes

    def get_animal_tissue_subgraph(self):
        tissues = self.get_animal_tissues()
        subgraph = self.get_graph().subgraph(tissues)
        return subgraph

if __name__ =='__main__':
    bto = BTO()
    bto_onto = networkx_ontology.NXOntology(bto.get_animal_tissue_subgraph())

    bto_onto.annotate_information_content()
    g = bto_onto.graph
    bodymap_tissues = [
    'BTO:0001487', 'BTO:0000047', 'BTO:0000142', 'BTO:0000149', 
    'BTO:0000269', 'BTO:0000562', 'BTO:0000671', 'BTO:0000751', 
    'BTO:0000759', 'BTO:0000763', 'BTO:0000784', 'BTO:0000975', 
    'BTO:0001129', 'BTO:0001103', 'BTO:0001363', 'BTO:0001379']
    
    import itertools
    for n0, n1 in itertools.combinations(bodymap_tissues, 2):
        print g.node[n0]['name'], g.node[n1]['name'], bto_onto.intrinsic_similarity(n0, n1)
    