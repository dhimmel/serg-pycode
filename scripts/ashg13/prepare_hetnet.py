import hetnet
import hetnet.agents
import hetnet.algorithms

network_dir = '/home/dhimmels/Documents/serg/ashg13/130904-1'
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
for cutoff in range(2, 5):
    identifier = 'length-{}-cutoff'.format(cutoff)
    metapaths_agent = hetnet.agents.MetaPathsAgent(metaedge_agent, identifier)
    metapaths = metagraph.extract_metapaths('gene', 'disease', cutoff)
    metapaths_agent.set(metapaths)
    metapaths = metapaths_agent.get()
    

# Create LearningEdges
positives = [edge for edge in graph.get_edges(False) if edge.metaedge == metaedge]
learning_edges = hetnet.algorithms.matched_negatives(positives)
learning_edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, 'matched-negative-1s-1t')
learning_edges_agent.set(learning_edges)
learning_edges = learning_edges_agent.get()

# Create Features
features_agent = hetnet.agents.FeaturesAgent(graph_agent, 'all-features')
features = hetnet.algorithms.get_features()
features_agent.set(features)
features = features_agent.get()

