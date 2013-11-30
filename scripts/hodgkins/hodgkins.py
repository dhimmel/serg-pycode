import itertools
import os
import csv

import networkx
import numpy
import numpy.linalg

import bioparser.data

import random_walk
import vegas_reader




def get_overlap(set_0, set_1):
    assert isinstance(set_0, set)
    assert isinstance(set_1, set)
    intersect = len(set_0 & set_1)
    union = len(set_0 | set_1)
    jaccard = float(intersect) / union if union else 0.0
    overlap = {'intersect': intersect, 'union': union, 'jaccard': jaccard}
    return overlap

gwas_catalog = bioparser.data.Data().gwas_catalog
doid_id_to_genes = gwas_catalog.get_doid_id_to_genes(p_cutoff=None, fdr_cutoff=None, mapped_term_cutoff=1, exclude_pmids=set())

doid = bioparser.data.Data().doid
doid.annotate_categories()
doid_ontology = bioparser.data.Data().doid.get_ontology()
#doid_ontology.graph

diseasome = networkx.Graph()

for doid_id, genes in doid_id_to_genes.items():
    data = doid_ontology.graph.node[doid_id].copy()
    data['genes'] = genes
    data['random_walk'] = dict()
    diseasome.add_node(doid_id, data)

for node_0, node_1 in itertools.combinations(diseasome, 2):
    data_0 = diseasome.node[node_0]
    data_1 = diseasome.node[node_1]
    genes_0 = data_0['genes']
    genes_1 = data_1['genes']
    edge_data = get_overlap(genes_0, genes_1)
    diseasome.add_edge(node_0, node_1, edge_data)


# Remove diseases without any gene overlap with other diseases
keep_nodes = list()
for node in diseasome:
    intersection_total = sum(data['intersect'] for data in diseasome[node].values())
    if intersection_total > 0:
        keep_nodes.append(node)

diseasome = diseasome.subgraph(keep_nodes)

hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'
vegas_dir = os.path.join(hodgkin_dir, 'vegas')
disease_to_genes = vegas_reader.genes_from_directory(vegas_dir, fdr_cutoff=0.8)

for hl_subtype, genes in disease_to_genes.iteritems():
    for node, data in diseasome.nodes(data=True):
        node_genes = data['genes']
        overlap = get_overlap(node_genes, genes)
        data['random_walk'][hl_subtype] = overlap

#gml_path = os.path.join(hodgkin_dir, 'diseasome-all.gml')
#networkx.write_gml(diseasome, gml_path)

jaccard_matrix = networkx.adjacency_matrix(diseasome, weight='jaccard')

resolvable_subtypes = list()
for hl_subtype in disease_to_genes.keys():
    jaccard_seed = numpy.array([data['random_walk'][hl_subtype]['jaccard'] for data in diseasome.node.values()])
    if not sum(jaccard_seed):
        print 'Skipping {} because all zero seeds'.format(hl_subtype)
        continue
    else:
        resolvable_subtypes.append(hl_subtype)
    jaccard_proximity, steps = random_walk.random_walk(r=0.2, seed_vector=jaccard_seed, adj_matrix=jaccard_matrix)
    jaccard_proximity = jaccard_proximity.T.tolist()[0]
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





