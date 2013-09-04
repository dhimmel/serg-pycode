import argparse
import collections
import csv
import itertools
import os

import hetnet
import hetnet.agents
import hetnet.algorithms


def analysis(identifier, network_dir, features_id, metaedge_id, metapaths_id, learning_edges_id):
    """ """
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    metaedge_agent = hetnet.agents.MetaEdgeAgent(graph_agent, metaedge_id)
    analysis_agent = hetnet.agents.AnalysisAgent(metaedge_agent, identifier)
    analysis_agent.set(features_id, metapaths_id, learning_edges_id)

    graph = graph_agent.get()
    features = analysis_agent.features_agent.get()
    metapaths = analysis_agent.metapaths_agent.get()
    learning_edges = analysis_agent.learning_edges_agent.get()
    
    result_gen = edge_feature_generator(graph, learning_edges, features, metapaths)
    
    feature_fieldnames = ['{}:{}'.format(feature['name'], metapath) 
        for feature, metapath in itertools.product(features, metapaths)]
    identity_fieldnames = ['source', 'target', 'status', 'group']
    fieldnames = identity_fieldnames + feature_fieldnames
    
    path = os.path.join(analysis_agent.directory, 'learning-features.txt')
    write_file = open(path, 'w')
    writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()
    
    for learning_edge, feature_dict in result_gen:
        source, target, kind, direction = learning_edge['edge_id']
        rowdict = collections.OrderedDict()
        rowdict['source'] = source
        rowdict['target'] = target
        rowdict['status'] = learning_edge['status']
        rowdict['group'] = learning_edge['group']
        rowdict.update(feature_dict)
        writer.writerow(rowdict)
        print source, target
        
    write_file.close()

def edge_feature_generator(graph, learning_edges, features, metapaths):
    """Returns a generator of (learning_edge, feature_dict) tuples"""
    for learning_edge in learning_edges:
        source, target, kind, direction = learning_edge['edge_id']
        source = graph.get_node(source)
        target = graph.get_node(target)
        exclude_edges = learning_edge['exclusions']
        feature_dict = features_betweens(graph, source, target, metapaths, features, exclude_edges)
        yield learning_edge, feature_dict

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
    

def features_betweens(graph, source, target, metapaths, features, exclude_edges=set()):
    results = collections.OrderedDict()
    for metapath in metapaths:
        feature_dict = features_between(graph, source, target, metapath, features, exclude_edges)
        for name, value in feature_dict.iteritems():
            key = '{}:{}'.format(name, metapath)
            results[key] = value
    return results

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
            paths['paths_st'] = [path for path in paths['paths_s'] if path.target() == target]
            
        elif 'paths_t' in path_types:
            paths_ts = [path for path in paths['paths_t'] if path.target() == source]
            paths['paths_st'] = [hetnet.Path(path.inverse_edges) for path in paths_ts]
            
        else:
            paths['paths_st'] = graph.paths_between(source, target, metapath, exclude_edges=exclude_edges)
        
    return paths


if __name__ == '__main__':
    network_dir = '/home/dhimmels/Documents/serg/ashg13/130904-1'
    features_id = 'all-features'
    metaedge_id = 'GaD-both'
    metapaths_id = 'length-2-cutoff'
    learning_edges_id = 'matched-negative-1s-1t'
    identifier = 'trial2'
    analysis(identifier, network_dir, features_id, metaedge_id, metapaths_id, learning_edges_id)
