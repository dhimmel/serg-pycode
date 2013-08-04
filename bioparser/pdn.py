import csv
import os

import networkx

import data


class PDN_networkx(object):
    
    def __init__(self, pdn_dir=None):
        if pdn_dir is None:
            pdn_dir = data.source_data_dir('pdn')
        self.pdn_dir = pdn_dir

    def sif_generator(self, path):
        with open(path) as sif_file:
            reader = csv.reader(sif_file, delimiter='\t')
            for row in reader:
                yield row
    
    def networkx_from_sif(self, path):
        graph = networkx.MultiGraph()
        for row in self.sif_generator(path):
            graph.add_edge(row[0], row[2], key=row[1])
        return graph

    def get_graph(self):
        if not hasattr(self, 'graph'):
            sif_path = os.path.join(self.pdn_dir, 'pdn.sif')
            self.graph = self.networkx_from_sif(sif_path)
        mapping = {node: 'DOID_' + node for node in self.graph.nodes()}
        networkx.relabel_nodes(self.graph, mapping, False)
        return self.graph


class PDN(object):
    
    def __init__(self, pdn_dir=None):
        """
        Paper and documentation:
        doi:10.1371/journal.pone.0022670
        
        Download cytoscape file from
        http://www3.nd.edu/~dial/plosone/diseasenetworks/
        
        Open cytoscape session, select all, export table of nodes as nodes.csv
        and table of edges as edges.csv.
        """
        if pdn_dir is None:
            pdn_dir = data.source_data_dir('pdn')
        self.pdn_dir = pdn_dir

    def read_table(self, name):
        path = os.path.join(self.pdn_dir, name + '.csv')
        with open(path) as f:
            reader = csv.DictReader(f)
            for row in reader:
                yield row

    def get_nodes(self):
        rows = self.read_table('nodes')
        id_to_name = dict()
        for row in rows:
            id = row['canonicalName']
            name = row['Column 2']
            id_to_name[id] = name
        return id_to_name
    
    def get_edges(self):
        rows = self.read_table('edges')
        for row in rows:
            id_0, id_1 = row['canonicalName'].split(' (pp) ')
            score = float(row['Column 3'])
            yield id_0, id_1, score
    
    def get_combordities(self):
        id_to_name = self.get_nodes()
        edges = self.get_edges()
        comorbid_tuples = list()
        for id_0, id_1, score in edges:
            names = sorted(id_to_name[id] for id in [id_0, id_1])
            comorbid_tuple = tuple(names + [score])
            comorbid_tuples.append(comorbid_tuple)
        comorbid_tuples.sort()
        return comorbid_tuples

    def write_comorbidities(self):
        path = os.path.join(self.pdn_dir, 'comorbidities.txt')
        with open(path, 'w') as wf:
            writer = csv.writer(wf, delimiter='\t')
            writer.writerows(self.get_combordities())


if __name__ == '__main__':
    pdn = PDN()
    pdn.write_comorbidities()
    