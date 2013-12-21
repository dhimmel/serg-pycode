import collections

import hetnet
import hetnet.agents
import hetnet.algorithms

network_dir = '/home/dhimmels/Documents/serg/ashg13/131219-1'
graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
print 'loading graph'
graph = graph_agent.get()
print 'graph loaded'
metagraph = graph.metagraph

#graph_agent.set(graph)


# MetaEdgeAgent
metaedge_agent = hetnet.agents.MetaEdgeAgent(graph_agent, 'GaD-both')
metaedge = graph.metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_agent.set(metaedge)
metaedge = metaedge_agent.get()


# MetaPathsAgent
for cutoff in range(2, 4):
    identifier = 'length-{}-cutoff'.format(cutoff)
    metapaths_agent = hetnet.agents.MetaPathsAgent(metaedge_agent, identifier)
    metapaths = metagraph.extract_metapaths('gene', 'disease', cutoff)
    metapaths_agent.set(metapaths)
    metapaths = metapaths_agent.get()

"""
# Create LearningEdges
metaedge_to_edges = graph.get_metaedge_to_edges()
positives = metaedge_to_edges[metaedge]
learning_edges = hetnet.algorithms.matched_negatives(positives)
learning_edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, 'matched-negative-1s-1t')
learning_edges_agent.set(learning_edges)
learning_edges = learning_edges_agent.get()
"""

# Find protein coding genes and diseases with at least one association
metanode_to_nodes = graph.get_metanode_to_nodes()
genes = metanode_to_nodes[metaedge.source]
diseases = metanode_to_nodes[metaedge.target]
genes = [gene for gene in genes if gene.data['locus_group'] == 'protein-coding gene']
diseases = [disease for disease in diseases if len(disease.get_edges(metaedge.inverse))]
print len(genes), 'genes'
print len(diseases), 'diseases'

"""
# create all_edges

all_edges = hetnet.algorithms.product_learning_edges(graph, metaedge, genes, diseases)
all_edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, 'all-edges')
all_edges_agent.set(all_edges)
all_edges = all_edges_agent.get()
"""

# Create Features
#features_agent = hetnet.agents.FeaturesAgent(graph_agent, 'all-features')
#features = hetnet.algorithms.get_features()
#features_agent.set(features)
#features = features_agent.get()


# Create a counter with associated number of genes for each disease
disease_counter = collections.Counter()
for disease in diseases:
     disease_counter[disease] = len(disease.get_edges(metaedge.inverse))

gene_minimum = 10
for disease, gene_number in disease_counter.items():
    if gene_number < gene_minimum:
        continue
    #disease_name = disease.data['name']
    edges = hetnet.algorithms.product_learning_edges(graph, metaedge, genes, [disease])
    edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, disease.id_)
    edges_agent.set(edges)

    analysis_agent = hetnet.agents.AnalysisAgent(metaedge_agent, disease.id_)
    analysis_agent.write_config(features_id='dwpc_0.5', metapaths_id='custom', learning_edges_id=disease.id_)


