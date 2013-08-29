import argparse
import collections

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
    
    features = hetnet.algorithms.get_features()

    

def PCs(paths):
    return len(paths['paths_s'])

def PCt(paths):
    return len(paths['paths_t'])

def PC(paths):
    return len(paths['paths_st'])

def NPC(paths):
    denom = len(paths['paths_s']) + len(paths['paths_t'])
    if denom:
        return 2.0 * len(paths['paths_st']) / denom
    else:
        return None

def DWPC(paths, damping_exponent):
    return hetnet.algorithms.degree_weighted_path_count(paths['paths_st'], damping_exponent)

algorithm_to_path_types = {
    'PCs': {'paths_s'},
    'PCt': {'paths_t'},
    'PC': {'paths_st'},
    'NPC': {'paths_st', 'paths_s', 'paths_t'},
    'DWPC': {'paths_st'}}

algorithm_to_function = {
    'PCs': PCs,
    'PCt': PCt,
    'PC': PC,
    'NPC': NPC,
    'DWPC': DWPC}
    



def features_between(graph, source, target, metapath, features, exclude_edges=set()):
    
    path_types = set()
    for feature in features:
        algorithm = feature['algorithm']
        path_types |= algorithm_to_path_types[algorithm]
    
    paths = paths_between(graph, source, target, metapath, path_types, exclude_edges)
    
    results = collections.OrderedDict()
    for feature in features:
        function = algorithm_to_function[feature['algorithm']]
        results[feature['name']] = function(paths, **feature['arguments'])
    return results

def paths_between(graph, source, target, metapath, path_types, exclude_edges):
    paths = dict()
    
    if 'paths_s' in path_types:
        paths['paths_s'] = graph.paths_from(source, metapath, exclude_edges=exclude_edges)
    
    if 'paths_t' in path_types:
        paths['paths_t'] = graph.paths_from(target, metapath.inverse, exclude_edges=exclude_edges)
    
    if 'paths_st' in path_types:
        
        if 'paths_s' in path_types:
            paths['paths_st'] = [path.target == target for path in paths['paths_s']]
            
        elif 'paths_t' in path_types:
            paths_ts = [path.target == source for path in paths['paths_t']]
            paths['paths_st'] = [hetnet.Path(path.inverse_edges) for path in paths_ts]
            
        else:
            paths['paths_st'] = graph.paths_between(source, target, metapath, exclude_edges=exclude_edges)
        
    return paths



