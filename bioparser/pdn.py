import csv
import os

import networkx

import data

def sif_generator(path):
    with open(path) as sif_file:
        reader = csv.reader(sif_file, delimiter='\t')
        for row in reader:
            yield row

def networkx_from_sif(path):
    graph = networkx.MultiGraph()
    for row in sif_generator(path):
        graph.add_edge(row[0], row[2], key=row[1])
    return graph

class PDN(object):
    
    def __init__(self, pdn_dir=None):
        if pdn_dir is None:
            pdn_dir = data.source_data_dir('pdn')
        self.pdn_dir = pdn_dir
    
    def get_graph(self):
        if not hasattr(self, 'graph'):
            sif_path = os.path.join(self.pdn_dir, 'pdn.sif')
            self.graph = networkx_from_sif(sif_path)
        mapping = {node: 'DOID_' + node for node in self.graph.nodes()}
        networkx.relabel_nodes(self.graph, mapping, False)
        return self.graph
        
if __name__ == '__main__':
    pdn = PDN()
    