import collections
import itertools
import os
import gzip
import csv

import hetnet
import hetnet.agents
import hetnet.algorithms


network_dir = '/home/dhimmels/Documents/serg/ashg13/140305-algosweep'
graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph



# Define Metrics
metrics = hetnet.algorithms.get_features()

# Define Metapaths
metapaths = metagraph.extract_metapaths('gene', 'disease', max_length=3)


# Define pairs
metaedge_GaD = metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metanode_to_nodes = graph.get_metanode_to_nodes()
genes = metanode_to_nodes[metaedge_GaD.source]
diseases = metanode_to_nodes[metaedge_GaD.target]
for node_list in genes, diseases:
    node_list.sort(key=lambda node: node.id_)

diseases_10plus = [disease for disease in diseases if len(disease.get_edges(metaedge_GaD.inverse)) >= 10]
disease_gene_pairs = list(itertools.product(diseases_10plus, genes))
total_edges = len(disease_gene_pairs)



feature_names = ['{}:{}'.format(metric['name'], metapath)
    for metric, metapath in itertools.product(metrics, metapaths)]
fieldnames = ['source', 'target', 'target_name', 'status'] + feature_names

feature_path = os.path.join(network_dir, 'features-edge-exclusions.txt.gz')
feature_file = gzip.open(feature_path, 'w')

writer = csv.DictWriter(feature_file, delimiter='\t', fieldnames=fieldnames)
writer.writeheader()

for i, (disease, gene) in enumerate(disease_gene_pairs):
    edge = graph.edge_dict.get((gene.id_, disease.id_, 'association', 'both'))
    status = 1 if edge else 0
    exclude_edges = {edge, edge.inverse} if status else set()
    results = hetnet.algorithms.features_betweens(graph, gene,
        disease, metapaths, metrics, exclude_edges=exclude_edges)
    results['source'] = gene
    results['target'] = disease
    results['target_name'] = disease.data['name']
    results['status'] = status

    percent = 100.0 * i / total_edges
    writer.writerow(results)
    print '{:.1f}% -  {:10}{}'.format(percent, gene, results['target_name'])

feature_file.close()
