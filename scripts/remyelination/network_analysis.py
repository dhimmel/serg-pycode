import collections
import os
import csv

import hetnet
import hetnet.algorithms
import hetnet.agents
import bioparser.data

import reader

network_dir = '/home/dhimmels/Documents/serg/remyelination/networks/140106'
graph_agent = hetnet.agents.GraphAgent(network_dir)
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph
metanode_to_nodes = graph.get_metanode_to_nodes()

target_metapath = metagraph.extract_metapaths('drug', 'protein', 1)[0]
target_metaedge = target_metapath[0]
drug_metanode = target_metapath.source()
protein_metanode = target_metapath.target()

drug_nodes = metanode_to_nodes[drug_metanode]
protein_nodes = metanode_to_nodes[protein_metanode]


features = list()

feature = collections.OrderedDict()
feature['name'] = 'PC'
feature['algorithm'] = 'PC'
feature['arguments'] = dict()
features.append(feature)

feature = collections.OrderedDict()
damping_exponent = 0.5
feature['name'] = 'DWPC_{}'.format(damping_exponent)
feature['algorithm'] = 'DWPC'
feature['arguments'] = {'damping_exponent': damping_exponent}
features.append(feature)

feature_fieldnames = list()
for protein in protein_nodes:
    for feature in [feature['name'] for feature in features]:
        feature_fieldnames.append('{}-{}'.format(feature, protein))


metaedge_to_edges = graph.get_metaedge_to_edges()
target_edges = metaedge_to_edges[target_metaedge]
sea_p_cutoffs = [10.0 ** -e for e in range(41)]

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
        for protein in protein_nodes:
            feature_dict = hetnet.algorithms.features_between(graph, drug, protein, target_metapath, features,
                duplicates=False, masked=False, exclude_nodes=set(), exclude_edges=set())
            fields.extend(feature_dict.values())
        writer.writerow(fields)
        print drug
    write_file.close()

print 'complete, unmasking now'
graph.unmask()