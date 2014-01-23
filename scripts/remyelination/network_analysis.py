import collections
import os
import csv
import gzip

import hetnet
import hetnet.algorithms
import hetnet.agents
import bioparser.data

import reader

network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140122-human'
graph_agent = hetnet.agents.GraphAgent(network_dir)
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph
metanode_to_nodes = graph.get_metanode_to_nodes()


metapath_protein = metagraph.extract_metapaths('drug', 'protein', 1)[0] # D-t-P
metapath_gene = metagraph.extract_metapaths('drug', 'gene', 3)[1] # D-t-P-p-G-i-G
metapath_c2 = metagraph.extract_metapaths('drug', 'c2.cp.all', 3)[0] # D-t-P-p-G-m-C2
metapath_c3 = metagraph.extract_metapaths('drug', 'c3.all', 3)[0] # D-t-P-p-G-m-C3
metapath_c5 = metagraph.extract_metapaths('drug', 'c5.all', 3)[0] # D-t-P-p-G-m-C5

drug_nodes = metanode_to_nodes[metapath_protein.source()]
protein_nodes = metanode_to_nodes[metapath_protein.target()]
gene_nodes = metanode_to_nodes[metapath_gene.target()]
gene_nodes = [gene_node for gene_node in gene_nodes if gene_node.data['locus_group'] == 'protein-coding gene']
c2_nodes = metanode_to_nodes[metapath_c2.target()]
c3_nodes = metanode_to_nodes[metapath_c3.target()]
c5_nodes = metanode_to_nodes[metapath_c5.target()]

# sort node lists
node_lists = [drug_nodes, protein_nodes, gene_nodes, c2_nodes, c3_nodes, c5_nodes]
for node_list in node_lists:
    node_list.sort(key=lambda l: l.id_)


# define metrics
metrics = list()

"""
feature = collections.OrderedDict()
feature['name'] = 'PC'
feature['algorithm'] = 'PC'
feature['arguments'] = dict()
features.append(feature)
"""
metric = collections.OrderedDict()
damping_exponent = 0.5
metric['name'] = 'DWPC_{}'.format(damping_exponent)
metric['algorithm'] = 'DWPC'
metric['arguments'] = {'damping_exponent': damping_exponent}
metrics.append(metric)
"""
feature = collections.OrderedDict()
feature['name'] = 'logP'
feature['algorithm'] = 'PWPC'
feature['arguments'] = {'probability_key': 'log_p_value'}
features.append(feature)
"""

feature_fieldnames = list()
metapath_info = [
    (metapath_protein, protein_nodes),
    (metapath_gene, gene_nodes),
    (metapath_c2, c2_nodes),
    (metapath_c3, c3_nodes),
    (metapath_c5, c5_nodes)]


feature_fieldnames = ['identifier', 'metric', 'metapath', 'target_id', 'target_name', 'target_description']
feature_rows = list()
for metapath, targets in metapath_info:
    for target in targets:
        for metric in [metric['name'] for metric in metrics]:
            feature_id = '{}|{}|{}'.format(metric, metapath, target)
            target_name = target.data.get('name')
            target_description = target.data.get('description')
            feature_row = feature_id, metric, metapath, target, target_name, target_description
            feature_rows.append(feature_row)

features_file = open(os.path.join(network_dir, 'features', 'features.txt'), 'w')
writer = csv.writer(features_file, delimiter='\t')
writer.writerow(feature_fieldnames)
writer.writerows(feature_rows)
features_file.close()

screen_reader = reader.ScreenReader()
screen_reader.binarize_screen()
smiles_to_compound = screen_reader.get_smiles_to_compound()
source_fieldnames = ['name', 'canonical_smiles', 'mean_MBP', 'SEM_MBP', 'mean_PDGFR', 'SEM_PDGFR', 'status']
compounds = list()
for drug in drug_nodes:
    compound = smiles_to_compound[drug.id_]
    compounds.append(compound)

compounds_file = open(os.path.join(network_dir, 'features', 'compounds.txt'), 'w')
writer = csv.DictWriter(compounds_file, fieldnames=source_fieldnames, delimiter='\t', extrasaction='ignore')
writer.writeheader()
writer.writerows(compounds)
compounds_file.close()



metaedge_to_edges = graph.get_metaedge_to_edges()
target_edges = metaedge_to_edges[metapath_protein[0]]
#sea_p_cutoffs = [10.0 ** -e for e in range(41)]
sea_p_cutoffs = [10e-5]

for sea_p_cutoff in sea_p_cutoffs:
    print 'SEA P-value cutoff:', sea_p_cutoff
    for edge in target_edges:
        if edge.data['p_value'] > sea_p_cutoff:
            edge.mask()

    path = os.path.join(network_dir, 'features', 'SEA-pval-cutoff-{:.0e}.txt.gz'.format(sea_p_cutoff))
    write_file = gzip.open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    fieldnames = ['drug'] + feature_fieldnames

    for drug in drug_nodes:
        fields = list()
        for metapath, targets in metapath_info:
            for target in targets:
                feature_dict = hetnet.algorithms.features_between(graph, drug, target, metapath, metrics,
                    duplicates=False, masked=False, exclude_nodes=set(), exclude_edges=set())
                fields.extend(feature_dict.values())
        writer.writerow(fields)
        print drug
    write_file.close()

print 'complete, unmasking now'
graph.unmask()