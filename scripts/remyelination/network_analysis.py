import collections
import os
import csv

import hetnet
import hetnet.algorithms
import hetnet.agents
import bioparser.data

import reader

network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140117'
graph_agent = hetnet.agents.GraphAgent(network_dir)
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph
metanode_to_nodes = graph.get_metanode_to_nodes()


metapath_protein = metagraph.extract_metapaths('drug', 'protein', 1)[0] # D-t-P
metapath_gene = metagraph.extract_metapaths('drug', 'gene', 3)[1] # D-t-P-p-G-i-G
metapath_c2 = metagraph.extract_metapaths('drug', 'c2.all', 3)[0] # D-t-P-p-G-m-C2
metapath_c5 = metagraph.extract_metapaths('drug', 'c5.all', 3)[0] # D-t-P-p-G-m-C5
metapath_c7 = metagraph.extract_metapaths('drug', 'c7.all', 3)[0] # D-t-P-p-G-m-C7

drug_nodes = metanode_to_nodes[metapath_protein.source()]
protein_nodes = metanode_to_nodes[metapath_protein.target()]
gene_nodes = metanode_to_nodes[metapath_gene.target()]
gene_nodes = [gene_node for gene_node in gene_nodes if gene_node.data['locus_group'] == 'protein-coding gene']
c2_nodes = metanode_to_nodes[metapath_c2.target()]
c5_nodes = metanode_to_nodes[metapath_c5.target()]
c7_nodes = metanode_to_nodes[metapath_c7.target()]

features = list()

"""
feature = collections.OrderedDict()
feature['name'] = 'PC'
feature['algorithm'] = 'PC'
feature['arguments'] = dict()
features.append(feature)
"""
feature = collections.OrderedDict()
damping_exponent = 0.5
feature['name'] = 'DWPC_{}'.format(damping_exponent)
feature['algorithm'] = 'DWPC'
feature['arguments'] = {'damping_exponent': damping_exponent}
features.append(feature)
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
    (metapath_c5, c5_nodes),
    (metapath_c7, c7_nodes)]

for metapath, targets in metapath_info:
    for target in targets:
        for feature in [feature['name'] for feature in features]:
            feature_fieldnames.append('{}|{}|{}'.format(feature, metapath, target))


metaedge_to_edges = graph.get_metaedge_to_edges()
target_edges = metaedge_to_edges[metapath_protein[0]]
sea_p_cutoffs = [10.0 ** -e for e in range(41)]
#sea_p_cutoffs = [1.0]

for sea_p_cutoff in sea_p_cutoffs:
    print 'SEA P-value cutoff:', sea_p_cutoff
    for edge in target_edges:
        if edge.data['p_value'] > sea_p_cutoff:
            edge.mask()

    path = os.path.join(network_dir, 'features', 'SEA-pval-cutoff-{:.0e}.txt'.format(sea_p_cutoff))
    write_file = open(path, 'w')
    writer = csv.writer(write_file, delimiter='\t')
    fieldnames = ['drug'] + feature_fieldnames
    writer.writerow(fieldnames)

    for drug in drug_nodes:
        fields = [drug]
        for metapath, targets in metapath_info:
            for target in targets:
                feature_dict = hetnet.algorithms.features_between(graph, drug, target, metapath, features,
                    duplicates=False, masked=False, exclude_nodes=set(), exclude_edges=set())
                fields.extend(feature_dict.values())
        writer.writerow(fields)
        print drug
    write_file.close()

print 'complete, unmasking now'
graph.unmask()