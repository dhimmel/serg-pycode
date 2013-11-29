import itertools
import os
import csv

import networkx
import numpy
import numpy.linalg

import bioparser.data

import random_walk



def read_vegas(path):
    symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    gene_to_vegasp = dict()
    for row in reader:
        vegas_gene = row['Gene']
        gene = symbol_to_gene.get(vegas_gene)
        if not gene:
            continue
        #assert gene not in gene_to_vegasp
        gene_to_vegasp[gene] = float(row['Pvalue'])
    read_file.close()
    return gene_to_vegasp


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

doid_ontology = bioparser.data.Data().doid.get_ontology()
#doid_ontology.graph

diseasome = networkx.Graph()

for doid_id, genes in doid_id_to_genes.items():
    doid_data = doid_ontology.graph.node[doid_id]
    data = {'name': doid_data['name']}
    data['genes'] = genes
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
vegas_path = os.path.join(hodgkin_dir, 'vegas', 'meta_USC_UC_IARC_updated_all.hg19_vegas_results')
gene_to_vegasp = read_vegas(vegas_path)
vegas_cutoff = 0.01
vegas_genes = set()
for gene, vegas_p in gene_to_vegasp.iteritems():
    if vegas_p <= vegas_cutoff:
        vegas_genes.add(gene)

for node, data in diseasome.nodes(data=True):
    node_genes = data['genes']
    overlap = get_overlap(node_genes, vegas_genes)
    data['overlap'] = overlap

gml_path = os.path.join(hodgkin_dir, 'diseasome-all.gml')
networkx.write_gml(diseasome, gml_path)

jaccard_seed = numpy.array([data['overlap']['jaccard'] for data in diseasome.node.values()])
jaccard_matrix = networkx.adjacency_matrix(diseasome, weight='jaccard')


jaccard_proximity, steps = random_walk.random_walk(r=0.2, seed_vector=jaccard_seed, adj_matrix=jaccard_matrix)
jaccard_proximity = jaccard_proximity.T.tolist()[0]
print 'total steps', steps

write_path = os.path.join(hodgkin_dir, 'proximity-all.txt')
write_file = open(write_path, 'w')
fieldnames = ['doid_code', 'name', 'jaccard_overlap', 'jaccard_rw_proximity']
writer = csv.DictWriter(write_file, fieldnames=fieldnames, delimiter='\t')
writer.writeheader()
for node, jaccard_rw_proximity in zip(diseasome, jaccard_proximity):
    data = diseasome.node[node]
    name = data['name']
    jaccard_overlap = data['overlap']['jaccard']
    row = {'doid_code': node, 'name': name,
           'jaccard_overlap': jaccard_overlap, 'jaccard_rw_proximity':jaccard_rw_proximity}
    writer.writerow(row)
write_file.close()





