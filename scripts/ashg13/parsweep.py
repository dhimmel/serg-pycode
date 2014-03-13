import collections
import itertools
import os
import gzip
import random

import hetnet
import hetnet.agents
import hetnet.algorithms
import utilities.floats


network_dir = '/home/dhimmels/Documents/serg/ashg13/140310-parsweep'
graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph


metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_DaG = metaedge_GaD.inverse
metaedge_GiG = metagraph.get_edge(('gene', 'gene', 'interaction', 'both'))
metaedge_GfG = metagraph.get_edge(('gene', 'gene', 'function', 'both'))
metaedge_DsD = metagraph.get_edge(('disease', 'disease', 'similarity', 'both'))
metaedge_GeT = metagraph.get_edge(('gene', 'tissue', 'expression', 'both'))
metaedge_TeG = metaedge_GeT.inverse
metaedge_TcD = metagraph.get_edge(('tissue', 'disease', 'cooccurrence', 'both'))
metaedge_DcT = metaedge_TcD.inverse


# (disease, disease, similarity, both)
metapath_GaDsD = metagraph.get_metapath((metaedge_GaD, metaedge_DsD))
metapath_GaDsDsD = metagraph.get_metapath((metaedge_GaD, metaedge_DsD, metaedge_DsD))
lin_thresholds = utilities.floats.sequence(0.0, 0.5, 0.01)

# (gene, gene, function, both) information
metapath_GfGaD = metagraph.get_metapath((metaedge_GfG, metaedge_GaD))
metapath_GfGfGaD = metagraph.get_metapath((metaedge_GfG, metaedge_GfG, metaedge_GaD))
probability_thresholds = utilities.floats.sequence(0.2, 1.0, 0.05)

# (gene, tissue, expression, both)
metapath_GeTeGaD = metagraph.get_metapath((metaedge_GeT, metaedge_TeG, metaedge_GaD))
metapath_GeTcD = metagraph.get_metapath((metaedge_GeT, metaedge_TcD))
log10_expr_thresholds = utilities.floats.sequence(0.5, 4.0, 0.1)

# Add (disease, tissue, cooccurrence, both)
metapath_GaDcTcD = metagraph.get_metapath((metaedge_GaD, metaedge_DcT, metaedge_TcD))
r_scaled_thresholds = utilities.floats.sequence(10, 45, 1)

metanode_to_nodes = graph.get_metanode_to_nodes()
genes = metanode_to_nodes[metaedge_GaD.source]
diseases = metanode_to_nodes[metaedge_GaD.target]
for node_list in genes, diseases:
    node_list.sort(key=lambda node: node.id_)

diseases_10plus = [disease for disease in diseases if len(disease.get_edges(metaedge_DaG)) >= 10]
disease_gene_pairs = list(itertools.product(diseases_10plus, genes))

dgs_tuples = [
    (disease, gene, int((gene.id_, disease.id_, 'association', 'both') in graph.edge_dict))
    for disease, gene in disease_gene_pairs]

random.seed(0)
dgs_tuples = [(d, g, s) for d, g, s in dgs_tuples if s or random.random() < 0.025]

def get_gte_mask(key, threshold):
    return lambda edge_data: edge_data[key] >= threshold

#DsD
features = list()
for lin_threshold in lin_thresholds:
    mask = get_gte_mask('lin_similarity', lin_threshold)
    feature = dict()
    feature['metapath'] = metapath_GaDsD
    feature['metaedge_to_mask'] = {metaedge_DsD: mask}
    feature['name'] = 'GaDsD_DsD={}'.format(lin_threshold)
    features.append(feature)

    feature = dict()
    feature['metapath'] = metapath_GaDsDsD
    feature['metaedge_to_mask'] = {metaedge_DsD: mask}
    feature['name'] = 'GaDsDsD_DsD={}'.format(lin_threshold)
    features.append(feature)


#GfG
for probability_threshold in probability_thresholds:
    mask = get_gte_mask('probability', probability_threshold)

    feature = dict()
    feature['metapath'] = metapath_GfGaD
    feature['metaedge_to_mask'] = {metaedge_GfG: mask}
    feature['name'] = 'GfGaD_GfG={}'.format(probability_threshold)
    features.append(feature)

    feature = dict()
    feature['metapath'] = metapath_GfGfGaD
    feature['metaedge_to_mask'] = {metaedge_GfG: mask}
    feature['name'] = 'GfGfGaD_GfG={}'.format(probability_threshold)
    features.append(feature)

#GeT
for log10_expr_threshold in log10_expr_thresholds:
    mask = get_gte_mask('log10_expr', log10_expr_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTeGaD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask, metaedge_TeG: mask}
    feature['name'] = 'GeTeGaD_GeT={}'.format(log10_expr_threshold)
    features.append(feature)


#DcT
for r_scaled_threshold in r_scaled_thresholds:
    mask = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GaDcTcD
    feature['metaedge_to_mask'] = {metaedge_DcT: mask, metaedge_TcD: mask}
    feature['name'] = 'GaDcTcD_DcT={}'.format(r_scaled_threshold)
    features.append(feature)

#GeT and TcD
for log10_expr_threshold, r_scaled_threshold in itertools.product(log10_expr_thresholds, r_scaled_thresholds):
    mask_GeT = get_gte_mask('log10_expr', log10_expr_threshold)
    mask_TcD = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTcD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask_GeT, metaedge_TcD: mask_TcD}
    feature['name'] = 'GeTcD_GeT={}_TcD={}'.format(log10_expr_threshold, r_scaled_threshold)
    features.append(feature)

feature_array = [[0 for i in xrange(len(features))] for j in xrange(len(dgs_tuples))]

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
        dwpc = hetnet.algorithms.DWPC({'paths_st': paths}, damping_exponent=0.5)
        feature_array[edge_index][feature_index] = dwpc

feature_path = os.path.join(network_dir, 'features-subset-fine.txt.gz')
feature_file = gzip.open(feature_path, 'w')
fieldnames = ['source', 'target', 'target_name', 'status'] + [f['name'] for f in features]
feature_file.write('\t'.join(fieldnames))
for (disease, gene, status), dwpc_list in zip(dgs_tuples, feature_array):
    feature_file.write('\n')
    row = [gene.id_, disease.id_, disease.data['name'], status] + dwpc_list
    feature_file.write('\t'.join(map(str, row)))

feature_file.close()
import collections
import itertools
import os
import gzip
import random

import hetnet
import hetnet.agents
import hetnet.algorithms
import utilities.floats


network_dir = '/home/dhimmels/Documents/serg/ashg13/140310-parsweep'
graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph


metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_DaG = metaedge_GaD.inverse
metaedge_GiG = metagraph.get_edge(('gene', 'gene', 'interaction', 'both'))
metaedge_GfG = metagraph.get_edge(('gene', 'gene', 'function', 'both'))
metaedge_DsD = metagraph.get_edge(('disease', 'disease', 'similarity', 'both'))
metaedge_GeT = metagraph.get_edge(('gene', 'tissue', 'expression', 'both'))
metaedge_TeG = metaedge_GeT.inverse
metaedge_TcD = metagraph.get_edge(('tissue', 'disease', 'cooccurrence', 'both'))
metaedge_DcT = metaedge_TcD.inverse


# (disease, disease, similarity, both)
metapath_GaDsD = metagraph.get_metapath((metaedge_GaD, metaedge_DsD))
metapath_GaDsDsD = metagraph.get_metapath((metaedge_GaD, metaedge_DsD, metaedge_DsD))
#lin_thresholds = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
lin_thresholds = utilities.floats.sequence(0.0, 0.5, 0.01)

# (gene, gene, function, both) information
metapath_GfGaD = metagraph.get_metapath((metaedge_GfG, metaedge_GaD))
metapath_GfGfGaD = metagraph.get_metapath((metaedge_GfG, metaedge_GfG, metaedge_GaD))
#probability_thresholds = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
probability_thresholds = utilities.floats.sequence(0.2, 1.0, 0.05)

# (gene, tissue, expression, both)
metapath_GeTeGaD = metagraph.get_metapath((metaedge_GeT, metaedge_TeG, metaedge_GaD))
metapath_GeTcD = metagraph.get_metapath((metaedge_GeT, metaedge_TcD))
#log10_expr_thresholds = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
log10_expr_thresholds = utilities.floats.sequence(0.5, 4.0, 0.1)

# Add (disease, tissue, cooccurrence, both)
metapath_GaDcTcD = metagraph.get_metapath((metaedge_GaD, metaedge_DcT, metaedge_TcD))
#r_scaled_thresholds = [15, 20, 25, 30, 35, 40, 45]
r_scaled_thresholds = utilities.floats.sequence(10, 45, 1)

metanode_to_nodes = graph.get_metanode_to_nodes()
genes = metanode_to_nodes[metaedge_GaD.source]
diseases = metanode_to_nodes[metaedge_GaD.target]
for node_list in genes, diseases:
    node_list.sort(key=lambda node: node.id_)

diseases_10plus = [disease for disease in diseases if len(disease.get_edges(metaedge_DaG)) >= 10]
disease_gene_pairs = list(itertools.product(diseases_10plus, genes))

dgs_tuples = [
    (disease, gene, int((gene.id_, disease.id_, 'association', 'both') in graph.edge_dict))
    for disease, gene in disease_gene_pairs]

random.seed(0)
dgs_tuples = [(d, g, s) for d, g, s in dgs_tuples if s or random.random() < 0.025]
#dgs_tuples = [(d, g, s) for d, g, s in dgs_tuples]

def get_gte_mask(key, threshold):
    return lambda edge_data: edge_data[key] >= threshold

#DsD
features = list()
for lin_threshold in lin_thresholds:
    mask = get_gte_mask('lin_similarity', lin_threshold)
    feature = dict()
    feature['metapath'] = metapath_GaDsD
    feature['metaedge_to_mask'] = {metaedge_DsD: mask}
    feature['name'] = 'GaDsD_DsD={}'.format(lin_threshold)
    features.append(feature)

    feature = dict()
    feature['metapath'] = metapath_GaDsDsD
    feature['metaedge_to_mask'] = {metaedge_DsD: mask}
    feature['name'] = 'GaDsDsD_DsD={}'.format(lin_threshold)
    features.append(feature)


#GfG
for probability_threshold in probability_thresholds:
    mask = get_gte_mask('probability', probability_threshold)

    feature = dict()
    feature['metapath'] = metapath_GfGaD
    feature['metaedge_to_mask'] = {metaedge_GfG: mask}
    feature['name'] = 'GfGaD_GfG={}'.format(probability_threshold)
    features.append(feature)

    feature = dict()
    feature['metapath'] = metapath_GfGfGaD
    feature['metaedge_to_mask'] = {metaedge_GfG: mask}
    feature['name'] = 'GfGfGaD_GfG={}'.format(probability_threshold)
    features.append(feature)

#GeT
for log10_expr_threshold in log10_expr_thresholds:
    mask = get_gte_mask('log10_expr', log10_expr_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTeGaD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask, metaedge_TeG: mask}
    feature['name'] = 'GeTeGaD_GeT={}'.format(log10_expr_threshold)
    features.append(feature)


#DcT
for r_scaled_threshold in r_scaled_thresholds:
    mask = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GaDcTcD
    feature['metaedge_to_mask'] = {metaedge_DcT: mask, metaedge_TcD: mask}
    feature['name'] = 'GaDcTcD_DcT={}'.format(r_scaled_threshold)
    features.append(feature)

#GeT and TcD
for log10_expr_threshold, r_scaled_threshold in itertools.product(log10_expr_thresholds, r_scaled_thresholds):
    mask_GeT = get_gte_mask('log10_expr', log10_expr_threshold)
    mask_TcD = get_gte_mask('r_scaled', r_scaled_threshold)
    feature = dict()
    feature['metapath'] = metapath_GeTcD
    feature['metaedge_to_mask'] = {metaedge_GeT: mask_GeT, metaedge_TcD: mask_TcD}
    feature['name'] = 'GeTcD_GeT={}_TcD={}'.format(log10_expr_threshold, r_scaled_threshold)
    features.append(feature)

feature_array = [[0 for i in xrange(len(features))] for j in xrange(len(dgs_tuples))]

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
        paths = graph.paths_between_tree(gene, disease, feature['metapath'],
            duplicates=False, masked=False, exclude_nodes={}, exclude_edges={})
        dwpc = hetnet.algorithms.DWPC({'paths_st': paths}, damping_exponent=0.5)
        feature_array[edge_index][feature_index] = dwpc

feature_path = os.path.join(network_dir, 'features-2.5percent.txt.gz')
feature_file = gzip.open(feature_path, 'w')
fieldnames = ['source', 'target', 'target_name', 'status'] + [f['name'] for f in features]
feature_file.write('\t'.join(fieldnames))
for (disease, gene, status), dwpc_list in zip(dgs_tuples, feature_array):
    feature_file.write('\n')
    row = [gene.id_, disease.id_, disease.data['name'], status] + dwpc_list
    feature_file.write('\t'.join(map(str, row)))

feature_file.close()
