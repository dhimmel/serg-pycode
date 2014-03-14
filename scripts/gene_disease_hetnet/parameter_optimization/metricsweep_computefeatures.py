import collections
import itertools
import os
import gzip
import csv
import random

import hetnet
import hetnet.agents
import hetnet.algorithms

project_dir = '/home/dhimmels/Documents/serg/gene-disease-hetnet'
network_dir = os.path.join(project_dir, 'networks', '140313-metricsweep')

graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph



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

diseases_10plus = [disease for disease in diseases if len(disease.get_edges(metaedge_GaD.inverse)) >= 10]
disease_gene_pairs = list(itertools.product(diseases_10plus, genes))
total_edges = len(disease_gene_pairs)



feature_names = ['{}:{}'.format(metric['name'], metapath)
    for metapath, metric in itertools.product(metapaths, metrics)]
fieldnames = ['source', 'target', 'target_name', 'status'] + feature_names

feature_path = os.path.join(network_dir, 'features.txt.gz')
feature_file = gzip.open(feature_path, 'w')

writer = csv.DictWriter(feature_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()

random.seed(0)
negative_prob = 0.025
print 'Negative Inclusion Probability: {}'.format(negative_prob)

for i, (disease, gene) in enumerate(disease_gene_pairs):
    edge = graph.edge_dict.get((gene.id_, disease.id_, 'association', 'both'))
    status = 1 if edge else 0
    if not status or random.random() > negative_prob:
        continue

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
