import itertools
import os
import csv

import networkx

import bioparser.data

import random_walk
import vegas_reader
from diseasome import DiseasomeMaker

hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'

def read_file(name, directory=hodgkin_dir):
    path = os.path.join(directory, name)
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        yield row
    read_file.close()



node_to_category = {row['doid_code']: row['category'] for row in read_file('doid-categories.txt')}

def write_proximity_table(path, graph, proximity_dict):
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    fieldnames = ['doid_code', 'name', 'category'] + proximity_dict.keys()
    writer.writerow(fieldnames)
    proximity_lists = zip(*proximity_dict.values())
    for (node, data), proximity_row in zip(graph.nodes_iter(data=True), proximity_lists):
        row = (node, data['name'], data['category']) + proximity_row
        writer.writerow(row)
    write_file.close()



def remove_one_proximities(diseasome, path):
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    fieldnames = ['source', 'source_name', 'source_category', 'target', 'target_name', 'target_category', 'proximity']
    writer.writerow(fieldnames)
    diseases = diseasome.nodes()
    for disease in diseases:
        disease_data = diseasome.node[disease]
        disease_category = disease_data['category']
        disease_name = disease_data['name']
        diseasome_minus = DiseasomeMaker.graph_minus(diseasome, [disease])
        unconnected_nodes = {node for node, degree in diseasome_minus.degree_iter() if not degree}
        diseasome_minus = DiseasomeMaker.graph_minus(diseasome_minus, unconnected_nodes)
        adjacency_matrix = DiseasomeMaker.get_adjacency_matrix(diseasome_minus)
        seed_vector = DiseasomeMaker.get_seed_distance(diseasome_minus, disease_data['genes'])
        rw_proximity, steps = random_walk.random_walk(r=0.2,
            seed_vector=seed_vector, adj_matrix=adjacency_matrix)
        rw_proximity = rw_proximity.tolist()

        for (target, target_data), proximity in zip(diseasome_minus.nodes(True), rw_proximity):
            row = (disease, disease_name, disease_category,
                   target, target_data['name'], target_data['category'], proximity)
            writer.writerow(row)
    write_file.close()





def multiple_seed_proximities(diseasome, disease_to_genes, path):
    disease_to_proximity = dict()
    for disease, genes in disease_to_genes.items():
        adjacency_matrix = DiseasomeMaker.get_adjacency_matrix(diseasome)
        seed_vector = DiseasomeMaker.get_seed_distance(diseasome, genes)
        if not sum(seed_vector):
            print 'Skipping {} because all zero seeds'.format(hltype)
            continue
        rw_proximity, steps = random_walk.random_walk(r=0.2,
            seed_vector=seed_vector, adj_matrix=adjacency_matrix)
        rw_proximity = rw_proximity.tolist()
        print disease, steps
        disease_to_proximity[disease] = rw_proximity

    write_proximity_table(path, diseasome, disease_to_proximity)




##################################################################################################
################# main ############################
gene_minimum = 10



symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
hodgkin_genes = {symbol_to_gene.get(row['gene']) for row in read_file('hodgkin-genes.txt')}
hodgkin_genes.discard(None)

doid_exclusions = {'DOID:0050589', # inflammatory bowel disease because of (UC and Crohn's)
                    #'DOID:557', # kidney disease because of DOID:784 (chronic kidney failure)
                    #'DOID:3620', # central nervous system cancer
                    'DOID:2914'} # immune system disease
doid_include = {'DOID:8567'} # Hodgkin's lymphoma

diseasome_maker = DiseasomeMaker()
diseasome_maker.set_defaults()
diseasome_maker.node_to_genes['DOID:8567'] |= hodgkin_genes # Add additional hodgkin genes
diseasome = diseasome_maker.get_graph(gene_minimum=gene_minimum, exclude=doid_exclusions, include=doid_include)
DiseasomeMaker.connect(diseasome)
DiseasomeMaker.add_node_attribute(diseasome, 'category', node_to_category)

## remove_one_proximities
path = os.path.join(hodgkin_dir, 'pairwise-proximities.txt')
remove_one_proximities(diseasome, path)


## multiple_seed_proximities for hodgkings
"""
vegas_dir = os.path.join(hodgkin_dir, 'vegas')
hltype_to_genes = vegas_reader.genes_from_directory(vegas_dir,
    method=vegas_reader.get_nominal_genes, cutoff=0.001)
hltype_to_genes['all']

path = os.path.join(hodgkin_dir, 'hlsubset-proximity-nominal-0.01.txt')
multiple_seed_proximities(diseasome, hltype_to_genes, path)
"""





#### Write diseasome with added HL-all
diseasome = diseasome_maker.get_graph(gene_minimum=gene_minimum, exclude=doid_exclusions, include=doid_include)
diseasome_maker.connect(diseasome)
diseasome_maker.add_node_attribute(diseasome, 'category', node_to_category, default='')
doid_code_to_name = dict()
for node, data in diseasome.nodes_iter(data=True):
    data['genes'] = ', '.join(gene.symbol for gene in data['genes'])
    data['doid_code'] = node
    doid_code_to_name[node] = data['name']
networkx.relabel_nodes(diseasome, doid_code_to_name, copy=False)
gml_path = os.path.join(hodgkin_dir, 'diseasome.gml')
networkx.write_gml(diseasome, gml_path)

