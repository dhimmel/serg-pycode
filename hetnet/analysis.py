import argparse
import collections
import csv
import itertools
import os

import hetnet
import hetnet.agents
import hetnet.algorithms


def analysis(analysis_agent):
    """ """
    features = analysis_agent.features_agent.get()
    metapaths = analysis_agent.metapaths_agent.get()
    learning_edges = list(analysis_agent.learning_edges_agent.get())
    graph = analysis_agent.graph_agent.get()
    
    
    feature_fieldnames = ['{}:{}'.format(feature['name'], metapath) 
        for feature, metapath in itertools.product(features, metapaths)]
    identity_fieldnames = ['source', 'target', 'status', 'group']
    fieldnames = identity_fieldnames + feature_fieldnames
    
    path = os.path.join(analysis_agent.directory, 'learning-features.txt')
    write_file = open(path, 'w')
    writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
    writer.writeheader()

    result_gen = edge_feature_generator(graph, learning_edges, features, metapaths)
    total_edges = len(learning_edges)
    for i, result in enumerate(result_gen):
        learning_edge, feature_dict = result
        source, target, kind, direction = learning_edge['edge_id']
        rowdict = collections.OrderedDict()
        rowdict['source'] = source
        rowdict['target'] = target
        rowdict['status'] = learning_edge['status']
        rowdict['group'] = learning_edge['group']
        rowdict.update(feature_dict)
        writer.writerow(rowdict)
        percent = 100.0 * i / total_edges
        print '{:.1f}% -  {:10}{}'.format(percent, source, target)
        
    write_file.close()

def edge_feature_generator(graph, learning_edges, features, metapaths):
    """Returns a generator of (learning_edge, feature_dict) tuples"""
    for learning_edge in learning_edges:
        source, target, kind, direction = learning_edge['edge_id']
        source = graph.get_node(source)
        target = graph.get_node(target)
        exclude_edges = learning_edge['exclusions']
        feature_dict = hetnet.algorithms.features_betweens(graph, source, target, metapaths, features, exclude_edges)
        yield learning_edge, feature_dict



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/ashg13/131219-1')
    parser.add_argument('--metaedge-id', default='GaD-both')
    parser.add_argument('--identifier', required=True)
    parser.add_argument('--config', action='store_true')
    parser.add_argument('--compute', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir

    graph_agent = hetnet.agents.GraphAgent(args.network_dir)
    metaedge_agent = hetnet.agents.MetaEdgeAgent(graph_agent, args.metaedge_id)
    analysis_agent = hetnet.agents.AnalysisAgent(metaedge_agent, args.identifier)

    if args.config:
        analysis_agent.write_config()
    
    if args.compute:
        analysis_agent.get()
        analysis(analysis_agent)

        

