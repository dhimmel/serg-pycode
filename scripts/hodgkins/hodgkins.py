import itertools
import os
import csv

import networkx
import numpy
import numpy.linalg

import bioparser.data

import random_walk
import vegas_reader





class DiseasomeLord(object):

    def __init__(self):
        pass

    def get_gene_annotations(self):
        # Load gwas_catalog and gene annotations.
        gwas_catalog = bioparser.data.Data().gwas_catalog
        node_to_genes = gwas_catalog.get_doid_id_to_genes(p_cutoff=None, fdr_cutoff=None, mapped_term_cutoff=1, exclude_pmids=set())
        return node_to_genes

    def get_ontology(self):
        # Load the disease ontology.
        doid = bioparser.data.Data().doid
        ontology = doid.get_ontology()
        return ontology

    def set_defaults(self):
        self.ontology = self.get_ontology()
        self.node_to_genes = self.get_gene_annotations()

    def get_diseasome(self, gene_minimum=10, exclude={}):
        """
        Returns nodes only
        gene_minimum - nodes with fewer than gene_minimum annotated genes are excluded
        exclude - set of nodes to exclude from the analysis.
        """
        diseasome = networkx.Graph()
        keep_data_keys = ['name']
        for node, genes in self.node_to_genes.items():
            if node in exclude:
                continue
            if len(genes) < gene_minimum:
                continue
            data = self.ontology.graph.node[node]
            data = {key: data[key] for key in keep_data_keys}
            data['genes'] = genes
            diseasome.add_node(node, data)

        return diseasome

    @staticmethod
    def get_overlap(set_0, set_1):
        assert isinstance(set_0, set)
        assert isinstance(set_1, set)
        intersect = len(set_0 & set_1)
        union = len(set_0 | set_1)
        jaccard = float(intersect) / union if union else 0.0
        overlap = {'intersect': intersect, 'union': union, 'jaccard': jaccard}
        return overlap

    @staticmethod
    def connect_diseasome(diseasome, weight_metric='jaccard'):
        for node_0, node_1 in itertools.combinations(diseasome, 2):
            data_0 = diseasome.node[node_0]
            data_1 = diseasome.node[node_1]
            genes_0 = data_0['genes']
            genes_1 = data_1['genes']
            edge_metrics = DiseasomeLord.get_overlap(genes_0, genes_1)
            weight = edge_metrics[weight_metric]
            if weight:
                diseasome.add_edge(node_0, node_1, weight=weight)

        # Remove diseases without any gene overlap with other diseases
        unconnected_nodes = {node for node, degree in diseasome.degree_iter() if not degree}
        diseasome.remove_nodes_from(unconnected_nodes)
        #diseasome = diseasome.subgraph(keep_nodes) # if want to leave original unmodified

        return diseasome

    @staticmethod
    def get_adjacency_matrix(diseasome):
        matrix = networkx.adjacency_matrix(diseasome)
        return matrix

    @staticmethod
    def get_seed_distance(diseasome, genes, weight_metric='jaccard'):
        genes = set(genes)
        seed_list = list()
        for data in diseasome.node.values():
            overlap = DiseasomeLord.get_overlap(genes, data['genes'])
            seed_list.append(overlap[weight_metric])
        seed = numpy.array(seed_list)
        return seed

diseasome_lord = DiseasomeLord()
diseasome_lord.set_defaults()



hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'
vegas_dir = os.path.join(hodgkin_dir, 'vegas')
hltype_to_genes = vegas_reader.genes_from_directory(vegas_dir, fdr_cutoff=0.5)

hl_exclusions = {'DOID:8567', # Hodgkin's lymphoma
              'DOID:0050589', # inflammatory bowel disease because of (UC and Crohn's)
              'DOID:557', # kidney disease because of DOID:784 (chronic kidney failure)
              'DOID:3620', # central nervous system cancer
              'DOID:2914'} # immune system disease

diseasome = diseasome_lord.get_diseasome(exclude=hl_exclusions)
diseasome_lord.connect_diseasome(diseasome)

valid_hltypes = list()
hltype_to_proximity = dict()
for hltype, genes in hltype_to_genes.items():
    adjacency_matrix = diseasome_lord.get_adjacency_matrix(diseasome)
    seed_vector = diseasome_lord.get_seed_distance(diseasome, genes)
    if not sum(seed_vector):
        print 'Skipping {} because all zero seeds'.format(hltype)
        continue
    valid_hltypes.append(hltype)
    rw_proximity, steps = random_walk.random_walk(r=0.2,
        seed_vector=seed_vector, adj_matrix=adjacency_matrix)
    rw_proximity = rw_proximity.tolist()
    print hltype, steps
    hltype_to_proximity[hltype][rw_proximity]






for node, data in diseasome.nodes_iter(data=True):
    del data['genes']
gml_path = os.path.join(hodgkin_dir, 'diseasome.gml')
networkx.write_gml(diseasome, gml_path)





"""
jaccard_matrix = networkx.adjacency_matrix(diseasome, weight='jaccard')

resolvable_subtypes = list()
for hl_subtype in hltype_to_genes.keys():
    jaccard_seed = numpy.array([data['random_walk'][hl_subtype]['jaccard'] for data in diseasome.node.values()])
    if not sum(jaccard_seed):
        print 'Skipping {} because all zero seeds'.format(hl_subtype)
        continue
    else:
        resolvable_subtypes.append(hl_subtype)
    jaccard_proximity, steps = random_walk.random_walk(r=0.2, seed_vector=jaccard_seed, adj_matrix=jaccard_matrix)
    for node, jaccard_rw_proximity in zip(diseasome, jaccard_proximity):
        diseasome.node[node]['random_walk'][hl_subtype]['jaccard_proximity'] = jaccard_rw_proximity
    print hl_subtype, 'total steps', steps


write_path = os.path.join(hodgkin_dir, 'proximity.txt')
write_file = open(write_path, 'w')

fieldnames = ['doid_code', 'name', 'categories'] + resolvable_subtypes
writer = csv.DictWriter(write_file, fieldnames=fieldnames, delimiter='\t')
writer.writeheader()
for node, jaccard_rw_proximity in zip(diseasome, jaccard_proximity):
    data = diseasome.node[node]
    name = data['name']
    categories = '; '.join(data['category_names'])
    row = {'doid_code': node, 'name': name, 'categories': categories}
    for subtype in resolvable_subtypes:
        row[subtype] = data['random_walk'][subtype]['jaccard_proximity']
    writer.writerow(row)
write_file.close()

"""

