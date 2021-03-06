import collections
import itertools
import os
import gzip
import random
import csv
import argparse

import hetnet
import hetnet.algorithms
import hetnet.readwrite
import utilities.floats


# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('--network-dir', type=os.path.expanduser, default=
    '~/Documents/serg/gene-disease-hetnet/networks/140614-thresholdsweep')
parser.add_argument('--partition-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
args = parser.parse_args()


print 'loading graph'
pkl_path = os.path.join(args.network_dir, 'graph', 'graph.pkl.gz')
graph = hetnet.readwrite.graph.read_pickle(pkl_path)
metagraph = graph.metagraph
print 'graph loaded'


metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_DaG = metaedge_GaD.inverse
metaedge_GeT = metagraph.get_edge(('gene', 'tissue', 'expression', 'both'))
metaedge_TeG = metaedge_GeT.inverse
metaedge_TlD = metagraph.get_edge(('tissue', 'disease', 'localization', 'both'))
metaedge_DlT = metaedge_TlD.inverse


# (gene, tissue, expression, both)
metapath_GeTeGaD = metagraph.get_metapath((metaedge_GeT, metaedge_TeG, metaedge_GaD))
metapath_GeTlD = metagraph.get_metapath((metaedge_GeT, metaedge_TlD))
log10_expr_thresholds = utilities.floats.sequence(0.5, 4.0, 0.1)

# Add (disease, tissue, localization, both)
metapath_GaDlTlD = metagraph.get_metapath((metaedge_GaD, metaedge_DlT, metaedge_TlD))
r_scaled_thresholds = utilities.floats.sequence(10, 45, 1)

metanode_to_nodes = graph.get_metanode_to_nodes()
genes = metanode_to_nodes[metaedge_GaD.source]
diseases = metanode_to_nodes[metaedge_GaD.target]
for node_list in genes, diseases:
    node_list.sort(key=lambda node: node.id_)


negative_prob = 0.005
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



def get_gte_mask(key, threshold):
    return lambda edge_data: edge_data[key] >= threshold


features = list()

#GeT
for log10_expr_threshold in log10_expr_thresholds:
    mask = get_gte_mask('log10_expr', log10_expr_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTeGaD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask, metaedge_TeG: mask}
    feature['name'] = 'GeTeGaD_GeT={}'.format(log10_expr_threshold)
    features.append(feature)

#DlT
for r_scaled_threshold in r_scaled_thresholds:
    mask = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GaDlTlD
    feature['metaedge_to_mask'] = {metaedge_DlT: mask, metaedge_TlD: mask}
    feature['name'] = 'GaDlTlD_DlT={}'.format(r_scaled_threshold)
    features.append(feature)

#GeT and TlD
for log10_expr_threshold, r_scaled_threshold in itertools.product(log10_expr_thresholds, r_scaled_thresholds):
    mask_GeT = get_gte_mask('log10_expr', log10_expr_threshold)
    mask_TlD = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTlD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask_GeT, metaedge_TlD: mask_TlD}
    feature['name'] = 'GeTlD_GeT={}_TlD={}'.format(log10_expr_threshold, r_scaled_threshold)
    features.append(feature)

feature_array = [[0 for i in xrange(len(features))] for j in xrange(len(dgs_tuples))]

damping_exponent = 0.4
metaedge_to_edges = graph.get_metaedge_to_edges()
for feature_index, feature in enumerate(features):
    graph.unmask()
    # mask
    for metaedge, mask in feature['metaedge_to_mask'].items():
        for edge in metaedge_to_edges[metaedge]:
            if not mask(edge.data):
                edge.mask()
    print feature['name']
    for edge_index, (disease, gene, status) in enumerate(dgs_tuples):
        if status:
            edge = graph.edge_dict[(gene.id_, disease.id_, 'association', 'both')]
            exclude_edges = {edge, edge.inverse}
        else:
            exclude_edges = set()
        paths = graph.paths_between_tree(gene, disease, feature['metapath'],
            duplicates=False, masked=False, exclude_nodes=set(), exclude_edges=exclude_edges)
        dwpc = hetnet.algorithms.DWPC(paths, damping_exponent=damping_exponent, exclude_edges=exclude_edges)
        feature_array[edge_index][feature_index] = dwpc

feature_path = os.path.join(args.network_dir,
    'features-exp{}.txt.gz'.format(damping_exponent))
feature_file = gzip.open(feature_path, 'w')
fieldnames = ['source', 'target', 'target_name', 'status'] + [f['name'] for f in features]
feature_file.write('\t'.join(fieldnames))
for (disease, gene, status), dwpc_list in zip(dgs_tuples, feature_array):
    feature_file.write('\n')
    row = [gene.id_, disease.id_, disease.data['name'], status] + dwpc_list
    feature_file.write('\t'.join(map(str, row)))

feature_file.close()
