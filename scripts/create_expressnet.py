import csv
import os

import bioparser.data
import heteronets.nxutils

from bioparser.metathesaurus import Concept


# Define schema for network
edge_tuples = [('drug', 'gene', 'upreg'),
               ('drug', 'gene', 'downreg'),
               ('disease', 'gene', 'upreg'),
               ('disease', 'gene', 'downreg'),
               ('drug', 'disease', 'indication'),
               ('gene', 'gene', 'function')]
kind_to_abbrev = {'drug': 'C', 'disease': 'D', 'gene': 'G',
                  'indication': 'i', 'upreg': 'u', 'downreg': 'd', 'function': 'f'}

network_id = '130227-1'
g = heteronets.nxutils.create_undirected_network(edge_tuples, kind_to_abbrev,
    name='expressnet', prepared=False, network_id=network_id)

dvd = bioparser.data.Data().dvd
###################
## Node Creation ##

# Create a disease node for each mapped dvd disease
diseases = set(dvd.get_disease_to_concept_id().values())
for disease in diseases:
    g.add_node(disease, kind='disease')

# Create a drug node for each mapped dvd drug
drugs = set(dvd.get_drug_to_concept_id().values())
for drug in drugs:
    g.add_node(drug, kind='drug')

# Add name as a data attribute to drug and disease nodes
meta = bioparser.data.Data().metathesaurus
with meta:
    concept_id_to_concept = meta.shelves['concepts']
    for node, data in g.nodes(True):
        if data['kind'] not in ['drug', 'disease']:
            continue
        data['name'] = concept_id_to_concept[node].name

# Create a gene node for every HGNC gene
hgnc = bioparser.data.Data().hgnc
genes = set(gene.symbol for gene in hgnc.get_genes())
for gene in genes:
    g.add_node(gene, kind='gene')

###################
## Edge Creation ##
direction_to_edge_kind = {'down': 'downreg', 'up': 'upreg'}

# Create drug-gene edges
for drug, gene, direction in dvd.get_drug_expression():
    g.add_edge(drug, gene, key=direction_to_edge_kind[direction])
    
# Create disease-gene edges
for disease, gene, direction in dvd.get_disease_expression():
    g.add_edge(disease, gene, key=direction_to_edge_kind[direction])

# Create indication edges
nci_code_to_concept_id = meta.get_source_code_to_concept_id('NCI')
mesh_id_to_concept_id = meta.get_source_code_to_concept_id('MSH')

indication_path = '/home/dhimmels/Documents/serg/ictnet/ictnet-creation/output-tables/tb_disease_drug_net.txt'
with open(indication_path) as f:
    reader = csv.reader(f, delimiter='\t')
    disease_drug_tuples = [tuple(row) for row in reader]
for mesh_drug, nci_disease in disease_drug_tuples:
    disease = nci_code_to_concept_id.get(nci_disease)
    drug = mesh_id_to_concept_id.get(mesh_drug)
    if not disease or not drug:
        continue
    if disease in diseases and drug in drugs:
        g.add_edge(drug, disease, key='indication')

# Remove unconnected nodes
heteronets.nxutils.remove_unconnected_nodes(g)

kind_to_edges = heteronets.nxutils.get_kind_to_edges(g)
indications = kind_to_edges[('drug', 'disease', 'indication')]
print '------------------Indications-------------------'
for node, neighbor, edge_kind in indications:
    print g.node[node]['name'], '\t', g.node[neighbor]['name']

###################
## Network Stats ##
print '------------------NetworkInfo-------------------'
heteronets.nxutils.print_node_kind_counts(g)
heteronets.nxutils.print_edge_kind_counts(g)


network_dir = os.path.join(args.ipadir, 'networks', args.network_id)

def save_as_pickle(g, network_dir):
    graph_dir = os.path.join(network_dir, 'graphs')
    if not os.path.exists(graph_dir):
        os.makedirs(graph_dir)
    pkl_path = os.path.join(graph_dir, 'graph.pkl')
    networkx.write_gpickle(g, pkl_path)
    print 'IPA network saved to', pkl_path

if not os.path.isdir(network_dir):
    os.mkdir(network_dir)
    
save_as_pickle(g, network_dir)
