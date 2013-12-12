import itertools
import os
import csv

import networkx

import bioparser.data

import random_walk
import vegas_reader
from diseasome import DiseasomeMaker

hodgkin_dir = '/home/dhimmels/Documents/serg/hodgkins/'

def read_file(name, directory=os.path.join(hodgkin_dir, 'input')):
    path = os.path.join(directory, name)
    read_file = open(path)
    reader = csv.DictReader(read_file, delimiter='\t')
    for row in reader:
        yield row
    read_file.close()


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

def get_hodgkin_genes():
    symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
    hodgkin_genes = {symbol_to_gene.get(row['gene']) for row in read_file('hodgkin-genes.txt')}
    hodgkin_genes.discard(None)
    return hodgkin_genes

def write_disease_to_genes(disease_to_genes, path):
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow(['disease', 'gene'])
    for disease, genes in disease_to_genes.items():
        genes_str = '|'.join(sorted(gene.symbol for gene in genes))
        writer.writerow([disease, genes_str])
    write_file.close()





##################################################################################################
################# main ############################
gene_minimum = 10

node_to_category = {row['doid_code']: row['category'] for row in read_file('doid-categories.txt')}
results_dir = os.path.join(hodgkin_dir, 'results')
if not os.path.exists(results_dir):
    os.mkdir(results_dir)


# Create diseasome with all diseases with over gene_minimum genes forcing HL inclusion
doid_exclusions = {'DOID:0050589', # inflammatory bowel disease because of (UC and Crohn's)
                    #'DOID:557', # kidney disease because of DOID:784 (chronic kidney failure)
                    #'DOID:3620', # central nervous system cancer
                    'DOID:2914'} # immune system disease
doid_include = {'DOID:8567'} # Hodgkin's lymphoma

doid_goto_doid = {'DOID:1612': 'DOID:3459', # breast cancer --> breast carcinoma
                  'DOID:9256': 'DOID:1520', # colorectal cancer --> colon carcinoma
                  'DOID:10283': 'DOID:10286' # prostate cancer --> prostate carcinoma
                 }

diseasome_maker = DiseasomeMaker()
diseasome_maker.set_gene_annotations(doid_goto_doid)


diseasome_maker.node_to_genes['DOID:8567'] |= get_hodgkin_genes() # Add additional hodgkin genes
diseasome = diseasome_maker.get_graph(gene_minimum=gene_minimum,
                                      exclude=doid_exclusions,
                                      include=doid_include)
DiseasomeMaker.connect(diseasome)
DiseasomeMaker.add_node_attribute(diseasome, 'category', node_to_category)
DiseasomeMaker.get_stats(diseasome)

# save network
gml_path = os.path.join(results_dir, 'diseasome.gml')
DiseasomeMaker.save_as_gml(diseasome, gml_path)
txt_path = os.path.join(results_dir, 'diseasome-flat.txt')
DiseasomeMaker.save_flat_txt(diseasome, txt_path)

## remove_one_proximities
path = os.path.join(results_dir, 'pairwise-proximities.txt')
remove_one_proximities(diseasome, path)


## multiple_seed_proximities for hodgkings
vegas_dir = os.path.join(hodgkin_dir, 'input', 'vegas')
hltype_to_genes = vegas_reader.genes_from_directory(vegas_dir,
    method=vegas_reader.get_nominal_genes, cutoff=0.001)
hltype_to_genes_path = os.path.join(hodgkin_dir, 'results', 'hl-subtype-vegas-genes.txt')
write_disease_to_genes(hltype_to_genes, hltype_to_genes_path)

diseasome_nohl = diseasome_maker.get_graph(gene_minimum=gene_minimum, exclude=doid_exclusions | {'DOID:8567'})
DiseasomeMaker.connect(diseasome_nohl)
DiseasomeMaker.add_node_attribute(diseasome_nohl, 'category', node_to_category)

path = os.path.join(results_dir, 'hlsubset-proximity.txt')
multiple_seed_proximities(diseasome_nohl, hltype_to_genes, path)

