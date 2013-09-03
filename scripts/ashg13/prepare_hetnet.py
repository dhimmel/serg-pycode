import hetnet
import hetnet.agents
import hetnet.algorithms

network_dir = '/home/dhimmels/Documents/serg/ashg13/130814-1'
graph_agent = hetnet.agents.GraphAgent(network_dir)

# Load graph
graph = graph_agent.get()
metagraph = graph.metagraph


# MetaEdgeAgent
metaedge_agent = hetnet.agents.MetaEdgeAgent(graph_agent, 'GaD-both')
metaedge = graph.metagraph.get_edge(('gene', 'disease', 'association', 'both'))
metaedge_agent.set(metaedge)
metaedge = metaedge_agent.get()


# MetaPathsAgent
metapaths_agent = hetnet.agents.MetaPathsAgent(metaedge_agent, 'length-4-cutoff')
metapaths = metagraph.extract_metapaths('gene', 'disease', 4)
metapaths_agent.set(metapaths)
metapaths = metapaths_agent.get()

# MetaPathsAgent
metapaths_agent = hetnet.agents.MetaPathsAgent(metaedge_agent, 'length-3-cutoff')
metapaths = metagraph.extract_metapaths('gene', 'disease', 3)
metapaths_agent.set(metapaths)
metapaths = metapaths_agent.get()


# Create LearningEdges
positives = [edge for edge in graph.get_edges(False) if edge.metaedge == metaedge]
learning_edges = hetnet.algorithms.matched_negatives(positives)
learning_edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, 'matched-negative-1s-1t')
learning_edges_agent.set(learning_edges)
learning_edges = learning_edges_agent.get()

# Create Features
features_agent = hetnet.agents.FeaturesAgent(network_dir, 'all-features')
features = hetnet.algorithms.get_features()
features_agent.set(features)
features = features_agent.get()

