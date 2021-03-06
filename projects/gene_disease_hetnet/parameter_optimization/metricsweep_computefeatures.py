import collections
import itertools
import os
import gzip
import csv
import argparse

import hetnet
import hetnet.readwrite
import hetnet.algorithms

# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--network-dir', type=os.path.expanduser, default=
    '~/Documents/serg/gene-disease-hetnet/networks/140614-metricsweep')
parser.add_argument('--partition-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
args = parser.parse_args()

print 'loading graph'
pkl_path = os.path.join(args.network_dir, 'graph', 'graph.pkl.gz')
graph = hetnet.readwrite.graph.read_pickle(pkl_path)
metagraph = graph.metagraph
print 'graph loaded'


# Define Metrics
metrics = hetnet.algorithms.get_metrics()

# Define Metapaths
metapaths = metagraph.extract_metapaths('gene', 'disease', max_length=3)


# Define pairs
metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metanode_to_nodes = graph.get_metanode_to_nodes()
genes = sorted(metanode_to_nodes[metaedge_GaD.source])
diseases = sorted(metanode_to_nodes[metaedge_GaD.target])
for node_list in genes, diseases:
    node_list.sort(key=lambda node: node.id_)

negative_prob = 0.05
print 'Negative Inclusion Probability: {}'.format(negative_prob)
dgs_tuples = list()
part_file = gzip.open(args.partition_path)
part_reader = csv.DictReader(part_file, delimiter='\t')
for part_row in part_reader:
    if part_row['part'] == 'test':
        continue
    status = part_row['status']
    percentile = float(part_row['percentile'])
    disease_node = graph.node_dict[part_row['disease_code']]
    gene_node = graph.node_dict[part_row['gene_symbol']]
    if part_row['status'] == 'HC_primary':
        dgs_tuple = disease_node, gene_node, 1
        dgs_tuples.append(dgs_tuple)
    if part_row['status'] == 'negative' and percentile <= negative_prob:
        dgs_tuple = disease_node, gene_node, 0
        dgs_tuples.append(dgs_tuple)
part_file.close()
total_edges = len(dgs_tuples)


feature_names = ['{}:{}'.format(metric['name'], metapath)
    for metapath, metric in itertools.product(metapaths, metrics)]
fieldnames = ['source', 'target', 'target_name', 'status'] + feature_names

feature_path = os.path.join(args.network_dir, 'features.txt.gz')
feature_file = gzip.open(feature_path, 'w')

writer = csv.DictWriter(feature_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()


for i, (disease, gene, status) in enumerate(dgs_tuples):
    edge = graph.edge_dict.get((gene.id_, disease.id_, 'association', 'both'))

    exclude_edges = {edge, edge.inverse} if status else set()

    results = collections.OrderedDict()
    results['source'] = gene
    results['target'] = disease
    results['target_name'] = disease.data['name']
    results['status'] = status

    for metapath in metapaths:
        paths_s = graph.paths_from(gene, metapath, exclude_edges=exclude_edges)
        paths_t = graph.paths_from(disease, metapath.inverse, exclude_edges=exclude_edges)
        paths = [path for path in paths_s if path.target() == disease]
        arg_dict = {'paths_s': paths_s, 'paths_t': paths_t, 'paths': paths,
                    'exclude_edges': exclude_edges}
        for metric in metrics:
            arguments = metric['arguments']
            for key in set(arg_dict) & set(arguments):
                arguments[key] = arg_dict[key]
            feature = metric['fxn'](**arguments)
            feature_key = '{}:{}'.format(metric['name'], metapath)
            results[feature_key] = feature

    percent = 100.0 * i / total_edges
    writer.writerow(results)
    print '{:.1f}% -  {:10}{}'.format(percent, gene, results['target_name'])

feature_file.close()
