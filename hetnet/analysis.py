import argparse

import hetnet
import hetnet.agents
import hetnet.algorithms


def analysis(network_dir, features_id, metaedge_id, metapaths_id, learning_edges_id):
    pass


graph_agent = hetnet.agents.GraphAgent(network_dir)
features_agent = hetnet.agents.FeaturesAgent(network_dir, features_id)
metaedge_agent = hetnet.agents.MetaEdgeAgent(graph_agent, metaedge_id)
metapaths_agent = hetnet.agents.MetaPathsAgent(metaedge_agent, metapaths_id)
learning_edges_agent = hetnet.agents.LearningEdgesAgent(metaedge_agent, learning_edges_id)


graph = graph_agent.get()
features = features_agent.get()
metapaths = metapaths_agent.get()
learning_edges = learning_edges_agent.get()

a

algorithm_to_paths = {
'PCs': {'paths_s'},
'PCt': {'paths_t'},
'PC': {'paths'},
'NPC': {'paths', 'paths_s', 'paths_t'},
'DWPC': {'paths'}}



