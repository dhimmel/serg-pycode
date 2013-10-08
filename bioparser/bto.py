import os
import csv

import networkx

import data
import networkx_ontology

class BTO(networkx_ontology.NXObo):
    
    def __init__(self, directory=None):
        if directory is None:
            directory = data.current_path('bto')
        obo_filename = 'BrendaTissueOBO.txt'
        json_filename = 'bto.networkx.json'
        pkl_filename = 'bto.networkx.pkl'
        self.keep_attributes = ['name', 'def_', 'synonym']
        super(BTO, self).__init__(directory, obo_filename, json_filename, pkl_filename)

    def get_animal_tissues(self):
        root = 'BTO:0001489' # whole body (descendant of animal)
        ontology = self.get_nx_ontology()
        return ontology.get_descendents(root)

    def get_animal_tissue_subgraph(self):
        tissues = self.get_animal_tissues()
        subgraph = self.get_graph().subgraph(tissues)
        return subgraph

if __name__ =='__main__':
    bto = BTO()
    print len(bto.get_animal_tissues()), 'animal tissues'

